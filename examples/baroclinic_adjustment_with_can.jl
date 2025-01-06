# # Baroclinic adjustment
#
# In this example, we simulate the evolution and equilibration of a baroclinically
# unstable front.
#
# ## Install dependencies
#
# First let's make sure we have all required packages installed.

# ```julia
# using Pkg
# pkg"add Oceananigans, ClimaOceanBiogeochemistry, CairoMakie"
# ```
using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients
using ClimaOceanBiogeochemistry.CarbonSystemSolvers.UniversalRobustCarbonSolver: UniversalRobustCarbonSystem
using ClimaOceanBiogeochemistry.CarbonSystemSolvers: CarbonSystemParameters, CarbonChemistryCoefficients

using Adapt
using KernelAbstractions: @kernel, @index
using Oceananigans.Utils: KernelParameters, launch!
using Printf
using CairoMakie

arch = CPU()

# ## Grid
# We use a three-dimensional channel that is periodic in the `x` direction:
Lx = 1000kilometers # east-west extent [m]
Ly = 1000kilometers # north-south extent [m]
Lz = 500    # depth [m]

Nx = 64
Ny = 64
Nz = 50

grid = RectilinearGrid(arch,
                       size     = (Nx, Ny, Nz),
                       x        = (0, Lx),
                       y        = (-Ly/2, Ly/2),
                       z        = (-Lz, 0),
                       topology = (
                            Periodic, Bounded, Bounded,
                        ))

# For air-sea CO₂ fluxes, we use a FluxBoundaryCondition for the "top" of the DIC tracer.
# We'll write a callback to calculate the flux every time step.
CO₂_flux = Field{Center, Center, Nothing}(grid)

dic_bcs  = FieldBoundaryConditions(
    top = FluxBoundaryCondition(CO₂_flux)
    )

## These are filled in compute_CO₂_flux!
ocean_CO₂ = Field{Center, Center, Nothing}(grid)
atmos_CO₂ = Field{Center, Center, Nothing}(grid)
pH        = Field{Center, Center, Nothing}(grid)
set!(pH, 8.0)

cmhr⁻¹_per_ms⁻¹ = 1 / 3.6e5 # conversion factor from cm/hr to m/s

## Supply some coefficients and external data for the CO₂ flux calculation
Base.@kwdef struct CO₂_flux_parameters
    surface_wind_speed   = 10. # ms⁻¹
    applied_pressure     = 0.0 # atm
    atmospheric_pCO₂     = 280e-6 # atm
    exchange_coefficient = 0.337 # cm hr⁻¹
    silicate_conc        = 15e-3 # mol m⁻³
    initial_pH_guess     = 8.0
    reference_density    = 1024.5 # kg m⁻³
    schmidt_dic_coeff0   = 2116.8
    schmidt_dic_coeff1   = 136.25 
    schmidt_dic_coeff2   = 4.7353 
    schmidt_dic_coeff3   = 9.2307e-2
    schmidt_dic_coeff4   = 7.555e-4
    schmidt_dic_coeff5   = 660.0
end
adapt_structure( 
    to, cp::CO₂_flux_parameters
    ) = CO₂_flux_parameters(
           adapt(to, cp.surface_wind_speed),
           adapt(to, cp.applied_pressure),
           adapt(to, cp.atmospheric_pCO₂),
           adapt(to, cp.exchange_coefficient),
           adapt(to, cp.silicate_conc),
           adapt(to, cp.initial_pH_guess),
           adapt(to, cp.reference_density),
           adapt(to, cp.schmidt_dic_coeff0),
           adapt(to, cp.schmidt_dic_coeff1),
           adapt(to, cp.schmidt_dic_coeff2),
	       adapt(to, cp.schmidt_dic_coeff3),
	       adapt(to, cp.schmidt_dic_coeff4),
           adapt(to, cp.schmidt_dic_coeff5),
)

"""
    compute_schmidt_dic(
        grid,
        schmidt_number_dic, 
        temperature, 
        schmidt_dic_coeff0,
        schmidt_dic_coeff1,
        schmidt_dic_coeff2,
        schmidt_dic_coeff3,
        schmidt_dic_coeff4,
        schmidt_dic_coeff5,
        )

Compute the Schmidt number for dissolved inorganic carbon (DIC) based on the temperature Θᶜ (in degrees Celsius).

# Arguments
- `grid::RectilinearGrid`: The model grid.
- `schmidt_number_dic::Field{Center, Center, Nothing}`: The computed Schmidt number for DIC.
- `temperature::Field{Center, Center, Nothing}`: Temperature in degrees Celsius.
- `schmidt_dic_coeff0::Float64`: Coefficient for the constant term (default: 2116.8).
- `schmidt_dic_coeff1::Float64`: Coefficient for the linear term (default: 136.25).
- `schmidt_dic_coeff2::Float64`: Coefficient for the quadratic term (default: 4.7353).
- `schmidt_dic_coeff3::Float64`: Coefficient for the cubic term (default: 9.2307e-2).
- `schmidt_dic_coeff4::Float64`: Coefficient for the quartic term (default: 7.555e-4).
- `schmidt_dic_coeff5::Float64`: Normalization coefficient (default: 660.0).
"""
@kernel function compute_schmidt_dic!(
    grid,
    schmidt_number_dic,
    temperature, 
    schmidt_dic_coeff0,
    schmidt_dic_coeff1,
    schmidt_dic_coeff2,
    schmidt_dic_coeff3,
    schmidt_dic_coeff4,
    schmidt_dic_coeff5,
    )

    i, j = @index(Global, NTuple)
    k = size(grid, 3)

    schmidt_number_dic[i, j, 1] =  ( 
         schmidt_dic_coeff0 - 
         schmidt_dic_coeff1 * temperature[i, j, k] + 
         schmidt_dic_coeff2 * temperature[i, j, k]^2 - 
         schmidt_dic_coeff3 * temperature[i, j, k]^3 + 
         schmidt_dic_coeff4 * temperature[i, j, k]^4 
    ) / schmidt_dic_coeff5
end

"""
    compute_piston_velocity(
        grid,
        piston_velocity,
        surface_wind_speed, 
        exchange_coefficient, 
        schmidt_number,
        )

Compute the piston velocity for gas exchange at the ocean surface.

# Arguments
- `grid::RectilinearGrid`: The model grid.
- `piston_velocity::Field{Center, Center, Nothing}`: The computed piston velocity.
- `surface_wind_speed::Field{Center, Center, Nothing}`: The wind speed at the ocean surface.
- `exchange_coefficient::Float64`: The gas exchange coefficient.
- `schmidt_number::Field{Center, Center, Nothing}`: The Schmidt number.
"""
@kernel function compute_piston_velocity!(
    grid,
    piston_velocity,
    surface_wind_speed,
    exchange_coefficient,
    schmidt_number
    )

    i, j = @index(Global, NTuple)
    
    piston_velocity[i, j, 1] = exchange_coefficient[i, j, 1] * 
                               surface_wind_speed[i, j, 1]^2 / 
                               sqrt(schmidt_number[i, j, 1])
end

"""
    solve_ocean_pCO₂!(
        grid,
        reference_density,
        ocean_pCO₂, 
        atmospheric_CO₂_solubility, 
        oceanic_CO₂_solubility,
        Θᶜ, Sᴬ, Δpᵦₐᵣ, Cᵀ, Aᵀ, Pᵀ, Siᵀ, pH, pCO₂ᵃᵗᵐ
        )

Compute the oceanic pCO₂ using the UniversalRobustCarbonSystem solver.

# Arguments
- `grid::RectilinearGrid`: The model grid.
- `reference_density::Float64`: The reference density of seawater.
- `ocean_pCO₂::Field{Center, Center, Nothing}`: The computed oceanic pCO₂.
- `atmospheric_CO₂_solubility::Field{Center, Center, Nothing}`: The solubility of CO₂ in the atmosphere.
- `oceanic_CO₂_solubility::Field{Center, Center, Nothing}`: The solubility of CO₂ in the ocean.
- `Θᶜ::Field{Center, Center, Nothing}`: Temperature in degrees Celsius.
- `Sᴬ::Field{Center, Center, Nothing}`: Salinity in PSU.
- `Δpᵦₐᵣ::Field{Center, Center, Nothing}`: Applied pressure in atm.
- `Cᵀ::Field{Center, Center, Nothing}`: Total dissolved inorganic carbon (DIC) in mol kg⁻¹.
- `Aᵀ::Field{Center, Center, Nothing}`: Total alkalinity (ALK) in mol kg⁻¹.
- `Pᵀ::Field{Center, Center, Nothing}`: Phosphate concentration in mol kg⁻¹.
- `Siᵀ::Field{Center, Center, Center}`: Silicate concentration in mol kg⁻¹.
- `pH::Field{Center, Center, Center}`: The computed pH.
- `pCO₂ᵃᵗᵐ::Field{Center, Center, Nothing}`: The partial pressure of CO₂ in the atmosphere.
"""
@kernel function solve_ocean_pCO₂!(
    grid,
    solver_params,
    reference_density,
    ocean_pCO₂, 
    atmospheric_CO₂_solubility, 
    oceanic_CO₂_solubility,
    temperature, 
    salinity, 
    applied_pressure_bar, 
    DIC, 
    ALK, 
    PO4, 
    pH, 
    atmosphere_pCO₂
    )
    
    i, j = @index(Global, NTuple)
    k = size(grid, 3)

    ## compute oceanic pCO₂ using the UniversalRobustCarbonSystem solver
    CarbonSolved = UniversalRobustCarbonSystem(
                        pH      = pH[i, j, 1], 
                        pCO₂ᵃᵗᵐ = atmosphere_pCO₂[i, j, 1],
                        Θᶜ      = temperature[i, j, k], 
                        Sᴬ      = salinity[i, j, k], 
                        Δpᵦₐᵣ   = applied_pressure_bar[i, j, 1], 
                        Cᵀ      = DIC[i, j, k]/reference_density, 
                        Aᵀ      = ALK[i, j, k]/reference_density, 
                        Pᵀ      = PO4[i, j, k]/reference_density, 
                        params  = solver_params,
    )

    ocean_pCO₂[i, j, 1]                 = CarbonSolved.pCO₂ᵒᶜᵉ
    atmospheric_CO₂_solubility[i, j, 1] = CarbonSolved.Pᵈⁱᶜₖₛₒₗₐ
    oceanic_CO₂_solubility[i, j, 1]     = CarbonSolved.Pᵈⁱᶜₖₛₒₗₒ
    pH[i, j, 1]                         = CarbonSolved.pH
end

"""
    compute_CO₂_flux(
        grid,
        CO₂_flux,
        piston_velocity,
        atmospheric_pCO₂,
        oceanic_pCO₂,
        atmospheric_CO₂_solubility,
        oceanic_CO₂_solubility,
        reference_density
    )
        
Compute the flux of CO₂ between the atmosphere and the ocean.

# Arguments
- `grid::RectilinearGrid`: The model grid.
- `reference_density::Float64`: The reference density of seawater.
- `CO₂_flux::Field{Center, Center, Nothing}`: The computed CO₂ flux.
- `piston_velocity::Field{Center, Center, Nothing}`: The piston velocity for gas exchange at the ocean surface.
- `atmospheric_pCO₂::Float64`: The partial pressure of CO₂ in the atmosphere.
- `oceanic_pCO₂::Field{Center, Center, Nothing}`: The partial pressure of CO₂ in the ocean.
- `atmospheric_CO₂_solubility::Field{Center, Center, Nothing}`: The solubility of CO₂ in the atmosphere.
- `oceanic_CO₂_solubility::Field{Center, Center, Nothing}`: The solubility of CO₂ in the ocean.

# Notes
The convention is that a positive flux is upwards (outgassing), and a negative flux is downwards (uptake).
"""
@kernel function compute_CO₂_flux!(
    grid,
    reference_density,
    CO₂_flux,
    piston_velocity,
    atmospheric_pCO₂,
    oceanic_pCO₂,
    atmospheric_CO₂_solubility,
    oceanic_CO₂_solubility,
    )
    i, j = @index(Global, NTuple)

    ## compute CO₂ flux (-ve for uptake, +ve for outgassing since convention is +ve upwards)
    CO₂_flux[i, j, 1] = - piston_velocity[i, j, 1] * (
                 atmospheric_pCO₂[i, j, 1] * atmospheric_CO₂_solubility[i, j, 1] - 
                 oceanic_pCO₂[i, j, 1]     * oceanic_CO₂_solubility[i, j, 1]
            ) * reference_density # Convert mol kg⁻¹ m s⁻¹ to mol m⁻² s⁻¹
end

# Use default carbon system parameters
solver_params = CarbonSystemParameters()

"""
    calculate_air_sea_carbon_exchange!(simulation)
    
Returns the tendency of DIC in the top layer due to air-sea CO₂ flux
using the piston velocity formulation of Wanninkhof (1992) and the
solubility/activity of CO₂ in seawater.
"""
@inline function calculate_air_sea_carbon_exchange!(simulation, solver_params)
    grid = simulation.model.grid

## get coefficients from CO₂_flux_parameters struct
## I really want the option to take these from the model
    (; surface_wind_speed, 
        applied_pressure,
        atmospheric_pCO₂, 
        exchange_coefficient, 
        silicate_conc, 
        initial_pH_guess, 
        reference_density,
        schmidt_dic_coeff0, 
        schmidt_dic_coeff1, 
        schmidt_dic_coeff2, 
        schmidt_dic_coeff3, 
        schmidt_dic_coeff4, 
        schmidt_dic_coeff5,
        ) = CO₂_flux_parameters()

    CO₂_flux = simulation.model.tracers.DIC.boundary_conditions.top.condition
    
    Siᵀ = Field{Center, Center, Center}(grid)
    set!(Siᵀ, silicate_conc)

    ## compute schmidt number for DIC
    schmidt_dic = Field{Center, Center, Nothing}(grid)
    set!(schmidt_dic, 0)

    temperature = simulation.model.tracers.T

    kernel_args = (
        grid,
        schmidt_dic, 
        temperature,
        schmidt_dic_coeff0,
        schmidt_dic_coeff1,
        schmidt_dic_coeff2,
        schmidt_dic_coeff3,
        schmidt_dic_coeff4,
        schmidt_dic_coeff5,
    )

    launch!(arch, 
	    grid, 
	    :xy, 
	    compute_schmidt_dic!, 
	    kernel_args...,
    )
    
    ## compute gas exchange coefficient/piston velocity and correct with Schmidt number
    wind_speed = Field{Center, Center, Nothing}(grid)
    set!(wind_speed, surface_wind_speed)

    average_exchange_coefficient = Field{Center, Center, Nothing}(grid)
    set!(average_exchange_coefficient, exchange_coefficient*cmhr⁻¹_per_ms⁻¹)

    piston_velocity = Field{Center, Center, Nothing}(grid)
    set!(piston_velocity, 0)

    kernel_args = (
        grid,
        piston_velocity,
        wind_speed,
        average_exchange_coefficient,
        schmidt_dic,
        )

    ## This is a 2d only kernel, so no need for custom kernel_size to handle slicing of 3d inputs
    launch!(arch, 
	    grid,
	    :xy,
	    compute_piston_velocity!, 
	    kernel_args...,
    )

    ## compute oceanic pCO₂ using the UniversalRobustCarbonSystem solver
    #atmos_CO₂ = Field{Center, Center, Nothing}(grid)
    set!(atmos_CO₂, atmospheric_pCO₂)

    applied_pressure_bar = Field{Center, Center, Nothing}(grid)
    set!(applied_pressure_bar, applied_pressure)

    salinity = simulation.model.tracers.S
    DIC = simulation.model.tracers.DIC
    ALK = simulation.model.tracers.ALK
    PO₄ = simulation.model.tracers.PO₄
    #oceanic_pCO₂               = Field{Center, Center, Nothing}(grid)
    atmospheric_CO₂_solubility = Field{Center, Center, Nothing}(grid)
    oceanic_CO₂_solubility     = Field{Center, Center, Nothing}(grid)
    
    kernel_args = (
        grid,
        solver_params,
        reference_density,
        ocean_CO₂, 
        atmospheric_CO₂_solubility, 
        oceanic_CO₂_solubility,
        temperature, 
        salinity, 
        applied_pressure_bar, 
        DIC, 
        ALK,
        PO₄, 
        pH, 
        atmos_CO₂,
        )

    launch!(arch,
	    grid,
	    :xy,
	    solve_ocean_pCO₂!,
	    kernel_args...,
    )

    ## compute CO₂ flux (-ve for uptake, +ve for outgassing since convention is +ve upwards)
    kernel_args = (
        grid,
        reference_density,
        CO₂_flux,
        piston_velocity,
        atmos_CO₂, 
        ocean_CO₂,
        atmospheric_CO₂_solubility, 
        oceanic_CO₂_solubility,
        )

    launch!(arch, 
	    grid, 
	    :xy, 
	    compute_CO₂_flux!, 
	    kernel_args...,
    )
    return nothing
end

# ## Model
# We built a `HydrostaticFreeSurfaceModel` with an `ImplicitFreeSurface` solver.
# Regarding Coriolis, we use a beta-plane centered at 45° South.
# We use the `SeawaterBuoyancy` model with a linear equation of state,
# that depends on temperature and salinity
model = HydrostaticFreeSurfaceModel(; grid,
            buoyancy = SeawaterBuoyancy(
                equation_of_state=LinearEquationOfState(
                    thermal_expansion = 2e-4,
                    haline_contraction = 8e-4),
                ),
            biogeochemistry = CarbonAlkalinityNutrients(; 
                grid,
                maximum_net_community_production_rate = 0.5/365.25day,
                ),
            closure = CATKEVerticalDiffusivity(),
            coriolis = BetaPlane(
                latitude = -45,
                ),
            tracers = (
                :T, :S, :e, :DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe
                ),
            boundary_conditions = (; 
                DIC=dic_bcs
                ),
            momentum_advection = WENO(),
            tracer_advection = WENO(),
            )

# We want to initialize our model with a baroclinically unstable front plus some small-amplitude
# noise.
"""
    ramp(y, Δy)

Linear ramp from 0 to 1 between -Δy/2 and +Δy/2.

For example:
```
            y < -Δy/2 => ramp = 0
    -Δy/2 < y < -Δy/2 => ramp = y / Δy
            y >  Δy/2 => ramp = 1
```
"""
ramp(y, Δy) = min(max(0, y/Δy + 1/2), 1)
nothing #hide

# We then use `ramp(y, Δy)` to construct an initial buoyancy configuration of a baroclinically
# unstable front. The front has a buoyancy jump `Δb` over a latitudinal width `Δy`.

N² = 4e-6 # [s⁻²] buoyancy frequency / stratification
M² = 8e-8 # [s⁻²] horizontal buoyancy gradient

Δy = 50kilometers # width of the region of the front
Δb = Δy * M²      # buoyancy jump associated with the front
ϵb = 1e-2 * Δb    # noise amplitude
α = model.buoyancy.model.equation_of_state.thermal_expansion
β = model.buoyancy.model.equation_of_state.haline_contraction
g = model.buoyancy.model.gravitational_acceleration

#bᵢ(x, y, z) = N² * z + Δb * ramp(y, Δy) + ϵb * randn()
T₀=10
Tᵢ(x, y, z) = T₀ + 1/(α*g) * 
                (N² * z + Δb * ramp(y, Δy) + 
                ϵb * randn())
S₀=34.5
Sᵢ(x, y, z) = (1 .- (Tᵢ(x, y, z) .- T₀)) .+ S₀

# Follow a similar method for BGC initial conditions
C₀ = 1.5 # Surface nutrient concentration
Cᵢ(x, y, z) = 1 .- (Tᵢ(x, y, z) .- T₀)/10 .+ C₀
A₀ = 1.8 # Surface nutrient concentration
Aᵢ(x, y, z) = 1 .- (Tᵢ(x, y, z) .- T₀)/10 .+ A₀
N₀ = 36 # Surface nutrient concentration
Nᵢ(x, y, z) = max(0, 
            ((1 .- (Tᵢ(x, y, z) .- T₀) * 3) .+ N₀) * 1e-3
            )
P₀ = 1.0 # Surface nutrient concentration
Pᵢ(x, y, z) = max(0,
            (1 .- (Tᵢ(x, y, z) .- T₀) .+ P₀) * 1e-3
            ) 
F₀ = 1.0 # Surface nutrient concentration
Fᵢ(x, y, z) = max(0,
            (1 .- (Tᵢ(x, y, z) .- T₀) .+ F₀) * 1e-6
            )
# Surface detritus concentration
DOP₀ = 0. 
 # Surface detritus concentration
POP₀ = 0.

# Plot the initial temperature and salinity profiles
set_theme!(Theme(fontsize = 12, linewidth=2))
x_points = grid.xᶜᵃᵃ[1:Nx,1,1]
y_points = grid.yᵃᶜᵃ[1:Ny,1,1]
z_points = grid.zᵃᵃᶜ[1:Nz,1,1]

fig = Figure(size = (1200, 1200))
axT = Axis(
    fig[1:2,1]; 
    xlabel = "y (m)", 
    ylabel = "Horizontal Temperature\nprofile [∘C]",
    )
axS = Axis(
    fig[3:4,1];
    xlabel = "y (m)", 
    ylabel = "Horizontal Salinity\nprofile [PSU]",
    )
axC = Axis(
    fig[5:6,1]; 
    xlabel = "y (m)", 
    ylabel = "Horizontal Carbon\nprofile [mol m⁻³]",
    )
axN = Axis(
    fig[7:8,1]; 
    xlabel = "y (m)", 
    ylabel = "Horizontal Nutrient\nprofile [mmol m⁻³]",
    )
for z = z_points[Nz:-10:1] 
    T_profile = [Tᵢ(x_points[1], y, z) for y in y_points]
    S_profile = [Sᵢ(x_points[1], y, z) for y in y_points]
    C_profile = [Cᵢ(x_points[1], y, z) for y in y_points]
    A_profile = [Aᵢ(x_points[1], y, z) for y in y_points]
    N_profile = [Nᵢ(x_points[1], y, z) for y in y_points]
    P_profile = [Pᵢ(x_points[1], y, z) for y in y_points]
    F_profile = [Fᵢ(x_points[1], y, z) for y in y_points]
    lines!(axT, y_points, T_profile[:,1], label="Temperature")
    lines!(axS, y_points, S_profile[:,1], label="Salinity")
    lines!(axC, y_points, C_profile[:,1], label="DIC")
    lines!(axC, y_points, A_profile[:,1], label="ALK")
    lines!(axN, y_points, N_profile[:,1]*1e3 , label="Nitrate (x1e-3)")
    lines!(axN, y_points, P_profile[:,1]*16e3, label="Phosphate (x16e-3)")
    lines!(axN, y_points, F_profile[:,1]*1e7 , label="Iron (x1e-7)")
end 
nothing #hide
fig

set!(model, 
    T   = Tᵢ,
    S   = Sᵢ,
    e   = 1e-6,
    DIC = Cᵢ,
    ALK = Aᵢ,
    PO₄ = Pᵢ,
    NO₃ = Nᵢ,
    DOP = DOP₀,
    POP = POP₀,
    Fe  = Fᵢ,
    )

# Now let's built a `Simulation`.
Δt₀ = 1minute
stop_time = 10days

simulation = Simulation(
    model, 
    Δt=Δt₀, 
    stop_time=stop_time,
    )

# We add a `TimeStepWizard` callback to adapt the simulation's time-step,
wizard = TimeStepWizard(
    cfl=0.2, 
    max_change=1.1, 
    max_Δt=30minutes,
    )
simulation.callbacks[:wizard] = Callback(
    wizard, 
    IterationInterval(20),
    )

# Also, we add a callback to print a message about how the simulation is going,
wall_clock = [time_ns()]

function print_progress(sim)
    @printf("[%05.2f%%] i: %d, t: %s, wall time: %s, max(u): (%6.3e, %6.3e, %6.3e) m/s, next Δt: %s\n",
            100 * (sim.model.clock.time / sim.stop_time),
            sim.model.clock.iteration,
            prettytime(sim.model.clock.time),
            prettytime(1e-9 * (time_ns() - wall_clock[1])),
            maximum(abs, sim.model.velocities.u),
            maximum(abs, sim.model.velocities.v),
            maximum(abs, sim.model.velocities.w),
            prettytime(sim.Δt))

    wall_clock[1] = time_ns()
    
    return nothing
end

simulation.callbacks[:print_progress] = Callback(
    print_progress, 
    IterationInterval(100),
    )

# Dont forget to add the callback to compute the CO₂ flux
add_callback!(simulation, calculate_air_sea_carbon_exchange!, parameters = solver_params)

# ## Diagnostics/Output

# Add some diagnostics. Here, we save the buoyancy, ``b``, at the edges of our domain as well as
# the zonal (``x``) average of buoyancy.

u, v, w = model.velocities

averages = NamedTuple(
    name => Average(
        model.tracers[name], dims=1
        ) for name in keys(model.tracers)
        )
filename = "baroclinic_adjustment_with_can"
save_fields_interval = 0.5day

slicers = (west = (1, :, :),
           east = (grid.Nx, :, :),
           south = (:, 1, :),
           north = (:, grid.Ny, :),
           bottom = (:, :, 1),
           top = (:, :, grid.Nz),
           )

for side in keys(slicers)
    indices = slicers[side]
    simulation.output_writers[side] = JLD2OutputWriter(
        model, model.tracers;
        filename = filename * "_$(side)_slice",
        schedule = TimeInterval(save_fields_interval),
        overwrite_existing = true,
        indices,
        )
end

outputs = (; CO₂_flux, ocean_CO₂, atmos_CO₂)
simulation.output_writers[:surface_bcs] = JLD2OutputWriter(
        model, outputs;
        filename = filename * "_surface_bcs_slice",
        schedule = TimeInterval(save_fields_interval),
        overwrite_existing = true,
        )

simulation.output_writers[:zonal] = JLD2OutputWriter(
        model, averages;
        filename = filename * "_zonal_average",
        schedule = TimeInterval(save_fields_interval),
        overwrite_existing = true,
    )

# Now let's run!

@info "Running the simulation..."

run!(simulation)

@info "Simulation completed in " * prettytime(
    simulation.run_wall_time,
    )

# ## Visualization

# Now we are ready to visualize our resutls! We use `CairoMakie` in this example.
# On a system with OpenGL `using GLMakie` is more convenient as figures will be
# displayed on the screen.

# We load the saved buoyancy output on the top, bottom, and east surface as `FieldTimeSeries`es.
sides = keys(slicers)
slice_filenames = NamedTuple(
    side => filename * "_$(side)_slice.jld2" for side in sides
    )

T_timeseries = (
    east   = FieldTimeSeries(slice_filenames.east,   "T"),
    north  = FieldTimeSeries(slice_filenames.north,  "T"),
    bottom = FieldTimeSeries(slice_filenames.bottom, "T"),
    top    = FieldTimeSeries(slice_filenames.top,    "T"),
    )

avg_T_timeseries = FieldTimeSeries(
    filename * "_zonal_average.jld2", "T",
    )

S_timeseries = (
    east   = FieldTimeSeries(slice_filenames.east,   "S"),
    north  = FieldTimeSeries(slice_filenames.north,  "S"),
    bottom = FieldTimeSeries(slice_filenames.bottom, "S"),
    top    = FieldTimeSeries(slice_filenames.top,    "S"),
    )

avg_S_timeseries = FieldTimeSeries(
    filename * "_zonal_average.jld2", "S",
    )

N_timeseries = (
    east   = FieldTimeSeries(slice_filenames.east,   "NO₃"),
    north  = FieldTimeSeries(slice_filenames.north,  "NO₃"),
    bottom = FieldTimeSeries(slice_filenames.bottom, "NO₃"),
    top    = FieldTimeSeries(slice_filenames.top,    "NO₃"),
    )

avg_N_timeseries = FieldTimeSeries(
    filename * "_zonal_average.jld2", "NO₃",
    )

P_timeseries = (
    east   = FieldTimeSeries(slice_filenames.east,   "PO₄"),
    north  = FieldTimeSeries(slice_filenames.north,  "PO₄"),
    bottom = FieldTimeSeries(slice_filenames.bottom, "PO₄"),
    top    = FieldTimeSeries(slice_filenames.top,    "PO₄"),
    )

avg_P_timeseries = FieldTimeSeries(
    filename * "_zonal_average.jld2", "PO₄",
    )

F_timeseries = (
    east   = FieldTimeSeries(slice_filenames.east,   "Fe"),
    north  = FieldTimeSeries(slice_filenames.north,  "Fe"),
    bottom = FieldTimeSeries(slice_filenames.bottom, "Fe"),
    top    = FieldTimeSeries(slice_filenames.top,    "Fe"),
    )

avg_F_timeseries = FieldTimeSeries(
    filename * "_zonal_average.jld2", "Fe",
    )

avg_pCO₂_timeseries = FieldTimeSeries(
    filename * "_surface_bcs_slice.jld2", "ocean_CO₂",
    )

times = avg_T_timeseries.times

climT = (8.983273690787108, 13.237066630286838) #1.1 .* extrema(avg_T_timeseries[:,:,:,1])
climS = (33.466822580627365, 36.513653646838634) #1.0 .* extrema(avg_S_timeseries[:,:,:,1])
climN = (0.030898256009215423, 0.0400419506843612) #1.0 .* extrema(avg_N_timeseries[:,:,:,1])
climP = (0.0, 0.0030144832060180573) #1.0 .* extrema(avg_P_timeseries[:,:,:,1])
climF = (0.0, 3.015480561290601e-6) #1.0 .* extrema(avg_F_timeseries[:,:,:,1])
climCO₂ = (180e-6, 380e-6)

kwargsT = (colorrange = climT, colormap = :balance)
kwargsS = (colorrange = climS, colormap = :haline)
kwargsN = (colorrange = climN, colormap = :Blues)
kwargsP = (colorrange = climP, colormap = :Purples)
kwargsF = (colorrange = climF, colormap = :Reds)
kwargspCO₂ = (colorrange = climCO₂, colormap = :PRGn)

nothing #hide

# We build the coordinates. We rescale horizontal coordinates so that they correspond to kilometers.

x, y, z = nodes(T_timeseries.east)

x = x .* 1e-3 # convert m -> km
y = y .* 1e-3 # convert m -> km

x_xz = repeat(x, 1, Nz)
y_xz_north = y[end] * ones(Nx, Nz)
z_xz = repeat(reshape(z, 1, Nz), Nx, 1)

x_yz_east = x[end] * ones(Ny, Nz)
y_yz = repeat(y, 1, Nz)
z_yz = repeat(reshape(z, 1, Nz), grid.Ny, 1)

x_xy = x
y_xy = y
z_xy_top = z[end] * ones(grid.Nx, grid.Ny)
z_xy_bottom = z[1] * ones(grid.Nx, grid.Ny)
nothing #hide

# Then we create a 3D axis. We use `zonal_slice_displacement` to control where the plot of the instantaneous
# zonal average flow is located.
set_theme!(Theme(fontsize = 18, linewidth=2))
fig3d = Figure(size = (2400, 1200))

zonal_slice_displacement = 1.2
axT = Axis3(fig3d[1:4, 1:4], aspect=(1, 1, 1/4),
           xlabel="x (km)", ylabel="y (km)", zlabel="z (m)",
           limits = (
            (x[1], zonal_slice_displacement * x[end]), 
            (y[1], y[end]), 
            (z[1], z[end]),
            ),
           elevation = 0.45, 
           azimuth = 6.8,
           xspinesvisible = false, 
           zgridvisible=false,
           protrusions=40,
           perspectiveness=0.7,
           )
axS = Axis3(fig3d[1:4, 6:9], aspect=(1, 1, 1/4),
           xlabel="x (km)", ylabel="y (km)", zlabel="z (m)",
           limits = (
            (x[1], zonal_slice_displacement * x[end]), 
            (y[1], y[end]), 
            (z[1], z[end]),
            ),
           elevation = 0.45, 
           azimuth = 6.8,
           xspinesvisible = false, 
           zgridvisible=false,
           protrusions=40,
           perspectiveness=0.7,
           )
axCO₂ = Axis3(fig3d[1:4, 11:14], aspect=(1, 1, 1/4),
           xlabel="x (km)", ylabel="y (km)", zlabel="z (m)",
           limits = (
            (x[1], zonal_slice_displacement * x[end]), 
            (y[1], y[end]), 
            (z[1], z[end]),
            ),
           elevation = 0.45, 
           azimuth = 6.8,
           xspinesvisible = false, 
           zgridvisible=false,
           protrusions=40,
           perspectiveness=0.7,
           )
axP = Axis3(fig3d[5:8, 1:4], aspect=(1, 1, 1/4),
           xlabel="x (km)", ylabel="y (km)", zlabel="z (m)",
           limits = (
            (x[1], zonal_slice_displacement * x[end]), 
            (y[1], y[end]), 
            (z[1], z[end]),
            ),
           elevation = 0.45, 
           azimuth = 6.8,
           xspinesvisible = false, 
           zgridvisible=false,
           protrusions=40,
           perspectiveness=0.7,
           )
axN = Axis3(fig3d[5:8, 6:9], aspect=(1, 1, 1/4),
           xlabel="x (km)", ylabel="y (km)", zlabel="z (m)",
           limits = (
            (x[1], zonal_slice_displacement * x[end]), 
            (y[1], y[end]), 
            (z[1], z[end]),
            ),
           elevation = 0.45, 
           azimuth = 6.8,
           xspinesvisible = false, 
           zgridvisible=false,
           protrusions=40,
           perspectiveness=0.7,
           )
axF = Axis3(fig3d[5:8, 11:14], aspect=(1, 1, 1/4),
           xlabel="x (km)", ylabel="y (km)", zlabel="z (m)",
           limits = (
            (x[1], zonal_slice_displacement * x[end]), 
            (y[1], y[end]), 
            (z[1], z[end]),
            ),
           elevation = 0.45, 
           azimuth = 6.8,
           xspinesvisible = false, 
           zgridvisible=false,
           protrusions=40,
           perspectiveness=0.7,
           )

nothing #hide

# We use Makie's `Observable` to animate the data. To dive into how `Observable`s work we
# refer to [Makie.jl's Documentation](https://makie.juliaplots.org/stable/documentation/nodes/index.html).

n = Observable(1)

# Now let's make a 3D plot of the buoyancy and in front of it we'll use the zonally-averaged output
# to plot the instantaneous zonal-average of the buoyancy.

T_slices = (east   = @lift(interior(
                T_timeseries.east[$n], 1, :, :),
                ),
            north  = @lift(interior(
                T_timeseries.north[$n], :, 1, :),
                ),
            bottom = @lift(interior(
                T_timeseries.bottom[$n], :, :, 1),
                ),
            top    = @lift(interior(
                T_timeseries.top[$n], :, :, 1),
                ))

avg_T = @lift interior(avg_T_timeseries[$n], 1, :, :)

S_slices = (east   = @lift(interior(
                S_timeseries.east[$n], 1, :, :),
                ),
            north  = @lift(interior(
                S_timeseries.north[$n], :, 1, :),
                ),
            bottom = @lift(interior(
                S_timeseries.bottom[$n], :, :, 1),
                ),
            top    = @lift(interior(
                S_timeseries.top[$n], :, :, 1),
                ))

avg_S = @lift interior(avg_S_timeseries[$n], 1, :, :)

pCO₂_slices = @lift interior(
                avg_pCO₂_timeseries[$n], :, :, 1, 1)

N_slices = (east   = @lift(interior(
                N_timeseries.east[$n], 1, :, :),
                ),
            north  = @lift(interior(
                N_timeseries.north[$n], :, 1, :),
                ),
            bottom = @lift(interior(
                N_timeseries.bottom[$n], :, :, 1),
                ),
            top    = @lift(interior(
                N_timeseries.top[$n], :, :, 1),
                ))

avg_N = @lift interior(avg_N_timeseries[$n], 1, :, :)

P_slices = (east   = @lift(interior(
                P_timeseries.east[$n], 1, :, :),
                ),
            north  = @lift(interior(
                P_timeseries.north[$n], :, 1, :),
                ),
            bottom = @lift(interior(
                P_timeseries.bottom[$n], :, :, 1),
                ),
            top    = @lift(interior(
                P_timeseries.top[$n], :, :, 1),
                ))

avg_P = @lift interior(avg_P_timeseries[$n], 1, :, :)

F_slices = (east   = @lift(interior(
                F_timeseries.east[$n], 1, :, :),
                ),
            north  = @lift(interior(
                F_timeseries.north[$n], :, 1, :),
                ),
            bottom = @lift(interior(
                F_timeseries.bottom[$n], :, :, 1),
                ),
            top    = @lift(interior(
                F_timeseries.top[$n], :, :, 1),
                ))

avg_F = @lift interior(avg_F_timeseries[$n], 1, :, :)

surface!(axT, x_yz_east, y_yz, z_yz;    
        color = T_slices.east, 
        kwargsT...,
        )
surface!(axT, x_xz, y_xz_north, z_xz;   
        color = T_slices.north, 
        kwargsT...,
        )
surface!(axT, x_xy, y_xy, z_xy_bottom ; 
        color = T_slices.bottom, 
        kwargsT...,
        )
surface!(axT, x_xy, y_xy, z_xy_top;     
        color = T_slices.top, 
        kwargsT...,
        )
        
surface!(axS, x_yz_east, y_yz, z_yz;    
        color = S_slices.east, 
        kwargsS...,
        )
surface!(axS, x_xz, y_xz_north, z_xz;   
        color = S_slices.north, 
        kwargsS...,
        )
surface!(axS, x_xy, y_xy, z_xy_bottom ; 
        color = S_slices.bottom, 
        kwargsS...,
        )
surface!(axS, x_xy, y_xy, z_xy_top;     
        color = S_slices.top, 
        kwargsS...,
        )

sfCO₂ = surface!(axCO₂, x_xy, y_xy, z_xy_top;     
        color = pCO₂_slices, 
        kwargspCO₂...,
        )

surface!(axP, x_yz_east, y_yz, z_yz;
        color = P_slices.east, 
        kwargsP...,
        )
surface!(axP, x_xz, y_xz_north, z_xz;
        color = P_slices.north, 
        kwargsP...,
        )
surface!(axP, x_xy, y_xy, z_xy_bottom ;
        color = P_slices.bottom, 
        kwargsP...,
        )
surface!(axP, x_xy, y_xy, z_xy_top;
        color = P_slices.top, 
        kwargsP...,
        )

surface!(axN, x_yz_east, y_yz, z_yz;
        color = N_slices.east, 
        kwargsN...,
        )
surface!(axN, x_xz, y_xz_north, z_xz;
        color = N_slices.north, 
        kwargsN...,
        )
surface!(axN, x_xy, y_xy, z_xy_bottom ;
        color = N_slices.bottom, 
        kwargsN...,
        )
surface!(axN, x_xy, y_xy, z_xy_top;
        color = N_slices.top, 
        kwargsN...,
        )

surface!(axF, x_yz_east, y_yz, z_yz;
        color = F_slices.east, 
        kwargsF...,
        )
surface!(axF, x_xz, y_xz_north, z_xz;
        color = F_slices.north, 
        kwargsF...,
        )
surface!(axF, x_xy, y_xy, z_xy_bottom ;
        color = F_slices.bottom, 
        kwargsF...,
        )
surface!(axF, x_xy, y_xy, z_xy_top;
        color = F_slices.top, 
        kwargsF...,
        )

sfT = surface!(axT, 
    zonal_slice_displacement .* x_yz_east, 
    y_yz, 
    z_yz; 
    color = avg_T, 
    kwargsT...,
    )

sfS = surface!(axS, 
    zonal_slice_displacement .* x_yz_east, 
    y_yz, 
    z_yz; 
    color = avg_S,
    kwargsS...,
    )

sfP = surface!(axP,
    zonal_slice_displacement .* x_yz_east, 
    y_yz, 
    z_yz; 
    color = avg_P,
    kwargsP...,
    )

sfN = surface!(axN,
    zonal_slice_displacement .* x_yz_east, 
    y_yz, 
    z_yz; 
    color = avg_N,
    kwargsN...,
    )

sfF = surface!(axF,
    zonal_slice_displacement .* x_yz_east, 
    y_yz, 
    z_yz; 
    color = avg_F,
    kwargsF...,
    )

contour!(axT, y, z, avg_T; 
         transformation = (:yz, zonal_slice_displacement * x[end]),
         levels = 15, 
         colorrange = climT,
         linewidth = 1, 
         color = :white,
         )

contour!(axS, y, z, avg_S; 
         transformation = (:yz, zonal_slice_displacement * x[end]),
         levels = 15, 
         colorrange = climS,
         linewidth = 1, 
         color = :white,
         )

contour!(axP, y, z, avg_P;
            transformation = (:yz, zonal_slice_displacement * x[end]),
            levels = 15, 
            colorrange = climP,
            linewidth = 1, 
            color = :white,
            )

contour!(axN, y, z, avg_N;
            transformation = (:yz, zonal_slice_displacement * x[end]),
            levels = 15, 
            colorrange = climN,
            linewidth = 1, 
            color = :white,
            )

contour!(axF, y, z, avg_F;
            transformation = (:yz, zonal_slice_displacement * x[end]),
            levels = 15, 
            colorrange = climF,
            linewidth = 1, 
            color = :white,
            )

Colorbar(
    fig3d[1:4, 5], 
    sfT, 
    label = "T [∘C]", 
    height = 200, 
    tellheight=false,
    )

Colorbar(
    fig3d[1:4, 10], 
    sfS, 
    label = "S [PSU]", 
    height = 200, 
    tellheight=false,
    )

Colorbar(    
    fig3d[1:4, 15], 
    sfCO₂, 
    ticks = (
        [180e-6, 230e-6, 280e-6, 330e-6, 380e-6],
        ["180","230","280","330","380"],
        ),
    label = "pCO₂ [uatm]", 
    height = 200, 
    tellheight=false,
    )

Colorbar(
    fig3d[5:8, 5], 
    sfP, 
    ticks = (
        [0, 1e-3, 2e-3, 3e-3],
        ["0","1","2","3"],
        ),
    label = "PO₄ [mmol m⁻³]", 
    height = 200, 
    tellheight=false,
    )

Colorbar(
    fig3d[5:8, 10], 
    sfN, 
    ticks = (
        [32e-3,34e-3,36e-3,38e-3,40e-3],
        ["32","34","36","38","40"],
        ),
    label = "NO₃ [mmol m⁻³]", 
    height = 200, 
    tellheight=false,
    )   

Colorbar(    
    fig3d[5:8, 15], 
    sfF, 
    ticks = (
        [0, 1e-6, 2e-6, 3e-6],
        ["0","1","2","3"],
        ),
    label = "Fe [nmol m⁻³]", 
    height = 200, 
    tellheight=false,
    )

title = @lift "Temperature and Salinity at t = " * string(round(times[$n] / day, digits=1)) * " days"

fig3d[1, 1:15] = Label(
    fig3d, 
    title; 
    fontsize = 24, 
    tellwidth = false, 
    padding = (0, 0, 0, 0),
    )

current_figure() # hide
fig3d

# Finally, we add a figure title with the time of the snapshot and then record a movie.

frames = 1:length(times)

CairoMakie.record(fig3d, filename * ".mp4", frames, framerate=8) do i
    msg = string("Plotting frame ", i, " of ", frames[end])
    print(msg * " \r")
    n[] = i
end
nothing #hide

# ![](baroclinic_adjustment.mp4)
