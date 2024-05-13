# # Illustration of using a carbon system solver
#
# This example illustrates how to use ClimaOceanBiogeochemistry's
# `UniversalRobustCarbonSystem` model in a 0-d context.

using ClimaOceanBiogeochemistry.CarbonSystemSolvers.UniversalRobustCarbonSolver: UniversalRobustCarbonSystem
using ClimaOceanBiogeochemistry.CarbonSystemSolvers: CarbonChemistryCoefficients

using Oceananigans
using Oceananigans.Units

using Printf
using CairoMakie

#
# We are going to simulate two bottles of soda, one opened and left in the fridge
#   the other opened and left on the counter to go flat.
#

# ## Model setup 
# The 0-d grid represents a 10cm bottle of soda
grid = RectilinearGrid(size = 1,
                       z = (-1meter/10, 0),
                       topology = (Flat, Flat, Bounded))

# For CO₂ exchange, we use a FluxBoundaryCondition for the "top" of the dissolved inorganic carbon (DIC) tracer.
# We'll write a callback later to calculate the flux every time step.
co₂_flux = Field{Center, Center, Nothing}(grid)
dic_bcs  = FieldBoundaryConditions(top = FluxBoundaryCondition(co₂_flux))

# Field arrays to store pCO₂ values filled in `compute_co₂_flux!`
soda_co₂ = Field{Center, Center, Nothing}(grid)
atmos_co₂ = Field{Center, Center, Nothing}(grid)

# Build function for CO₂ flux calculation. Dissolved CO₂ in the soda will exchange with the overlying atmosphere
# These are some coefficients and constants that we'll use in the calculation
Base.@kwdef struct Pᶠˡᵘˣᶜᵒ²{FT}
    surface_wind_speed   :: FT = 10. # ms⁻¹
    atmospheric_pCO₂     :: FT = 280e-6 # atm
    exchange_coefficient :: FT = 0.337 # cm hr⁻¹
    salinity             :: FT = 35.0 # psu
    alkalinity           :: FT = 2.35e-3 # mol kg⁻¹
    silicate_conc        :: FT = 15e-6 # mol kg⁻¹
    phosphate_conc       :: FT = 0e-6 # mol kg⁻¹
    initial_pH_guess     :: FT = 8.0
    reference_density    :: FT = 1024.5 # kg m⁻³
    schmidt_dic_coeff0   :: FT = 2116.8
    schmidt_dic_coeff1   :: FT = 136.25 
    schmidt_dic_coeff2   :: FT = 4.7353 
    schmidt_dic_coeff3   :: FT = 9.2307e-2
    schmidt_dic_coeff4   :: FT = 7.555e-4
    schmidt_dic_coeff5   :: FT = 660.0
end

"""
    compute_co₂_flux!(simulation)
    
Returns the tendency due to CO₂ flux using the piston velocity 
formulation of Wanninkhof (1992) and the solubility/activity of 
CO₂ that depends on temperature, etc.
"""
@inline function compute_co₂_flux!(simulation)
## Get coefficients from Pᶠˡᵘˣᶜᵒ² struct
## I really want the option to take these from the model
    (; surface_wind_speed, 
       atmospheric_pCO₂, 
       exchange_coefficient, 
       salinity,
       alkalinity,
       silicate_conc, 
       phosphate_conc,
       initial_pH_guess, 
       reference_density,   
       schmidt_dic_coeff0, 
       schmidt_dic_coeff1, 
       schmidt_dic_coeff2, 
       schmidt_dic_coeff3, 
       schmidt_dic_coeff4, 
       schmidt_dic_coeff5,
       ) = Pᶠˡᵘˣᶜᵒ²()

    U₁₀      = surface_wind_speed  
    pCO₂ᵃᵗᵐ  = atmospheric_pCO₂    
    Kʷₐᵥₑ    = exchange_coefficient
    Sᴬ       = salinity
    Aᵀ       = alkalinity
    Siᵀ      = silicate_conc   
    Pᵀ       = phosphate_conc   
    pH       = initial_pH_guess
    ρʳᵉᶠ     = reference_density 
    a₀ˢᶜᵈⁱᶜ   = schmidt_dic_coeff0 
    a₁ˢᶜᵈⁱᶜ   = schmidt_dic_coeff1 
    a₂ˢᶜᵈⁱᶜ   = schmidt_dic_coeff2 
    a₃ˢᶜᵈⁱᶜ   = schmidt_dic_coeff3 
    a₄ˢᶜᵈⁱᶜ   = schmidt_dic_coeff4 
    a₅ˢᶜᵈⁱᶜ   = schmidt_dic_coeff5 
    co₂_flux = simulation.model.tracers.DIC.boundary_conditions.top.condition

    cmhr⁻¹_per_ms⁻¹ = 1 / 3.6e5 # conversion factor from cm/hr to m/s

    ## applied pressure in atm (i.e. Pˢᵘʳᶠ-Pᵃᵗᵐ) 
    ## Positive when the can is sealed, then zero when the can is opens
    ## On average, the 12 ounce soda cans sold in the US tend to have a pressure of roughly 120 kPa when canned at 4 °C
    if iteration(simulation) <= 1
        Δpᵦₐᵣ   = 0.2
    else
        ## *pssshhhht* the bottle is opened
        Δpᵦₐᵣ   = 0.0
    end

    ## access model fields
    Θᶜ = simulation.model.tracers.T[
        simulation.model.grid.Nx,
        simulation.model.grid.Ny,
        simulation.model.grid.Nz
        ]
    Cᵀ = simulation.model.tracers.DIC[
        simulation.model.grid.Nx,
        simulation.model.grid.Ny,
        simulation.model.grid.Nz
        ]/ρʳᵉᶠ # Convert mol m⁻³ to mol kg⁻¹

    ## compute soda pCO₂ using the UniversalRobustCarbonSystem solver
    ## Returns soda pCO₂ (in atm) and atmosphere/soda solubility coefficients (mol kg⁻¹ atm⁻¹)
    (; pCO₂ᵒᶜᵉ, Pᵈⁱᶜₖₛₒₗₐ, Pᵈⁱᶜₖₛₒₗₒ) = 
        UniversalRobustCarbonSystem(
                Θᶜ, Sᴬ, Δpᵦₐᵣ, Cᵀ, Aᵀ, Pᵀ, Siᵀ, pH, pCO₂ᵃᵗᵐ,
                )

    ## store the soda and atmospheric CO₂ concentrations into Fields
    soda_co₂[
        simulation.model.grid.Nx,
        simulation.model.grid.Ny,
        simulation.model.grid.Nz,
        ] = (pCO₂ᵒᶜᵉ * Pᵈⁱᶜₖₛₒₗₒ ) * ρʳᵉᶠ # Convert mol kg⁻¹ m s⁻¹ to mol m⁻² s⁻¹
    atmos_co₂[
        simulation.model.grid.Nx,
        simulation.model.grid.Ny,
        simulation.model.grid.Nz,
        ] = (pCO₂ᵃᵗᵐ * Pᵈⁱᶜₖₛₒₗₐ) * ρʳᵉᶠ # Convert mol kg⁻¹ to mol m⁻³ 

    ## compute schmidt number for DIC
    Cˢᶜᵈⁱᶜ =  ( a₀ˢᶜᵈⁱᶜ - 
                a₁ˢᶜᵈⁱᶜ * Θᶜ + 
                a₂ˢᶜᵈⁱᶜ * Θᶜ^2 - 
                a₃ˢᶜᵈⁱᶜ * Θᶜ^3 + 
                a₄ˢᶜᵈⁱᶜ * Θᶜ^4 
              )/a₅ˢᶜᵈⁱᶜ
    
    ## compute gas exchange coefficient/piston velocity and correct with Schmidt number
    Kʷ =  (
           (Kʷₐᵥₑ * cmhr⁻¹_per_ms⁻¹) * U₁₀^2
          ) / sqrt(Cˢᶜᵈⁱᶜ)    
    
    ## compute co₂ flux (-ve for uptake, +ve for outgassing since convention is +ve upwards in the soda)
    co₂_flux[
        simulation.model.grid.Nx,
        simulation.model.grid.Ny,
        simulation.model.grid.Nz
        ] = - Kʷ * (
                    pCO₂ᵃᵗᵐ * Pᵈⁱᶜₖₛₒₗₐ - 
                    pCO₂ᵒᶜᵉ * Pᵈⁱᶜₖₛₒₗₒ
                   ) * ρʳᵉᶠ # Convert mol kg⁻¹ m s⁻¹ to mol m⁻² s⁻¹
    return nothing
end

# # Simulation of a soda outgassing CO₂ in the fridge
model_open_in_the_fridge = NonhydrostaticModel(; 
    grid,
    velocities = nothing,
    buoyancy   = nothing,
    closure    = nothing,
    tracers    = (:T, :DIC),
    boundary_conditions = (; DIC=dic_bcs),
    )

# Initial conditions for the refridgerated soda
Tᵢ   = 4     # °C
DICᵢ = 2.4   # mol/m³ 
set!(model_open_in_the_fridge, T = Tᵢ, DIC = DICᵢ)

simulation = Simulation(model_open_in_the_fridge, Δt=10minutes, stop_time=24hours)

# Add an output writer...
fnm_fridge = "soda_in_the_fridge.jld2"
outputs = (; co₂_flux, 
             soda_co₂, 
             atmos_co₂, 
             model_open_in_the_fridge.tracers.T, 
             model_open_in_the_fridge.tracers.DIC,
             )
simulation.output_writers[:jld2] = JLD2OutputWriter(
    model_open_in_the_fridge, outputs;
    filename = fnm_fridge,
    schedule = TimeInterval(10minutes),
    overwrite_existing = true,
    )

# ... and don't forget to add a callback to compute the CO₂ flux
add_callback!(simulation, compute_co₂_flux!)

# Run the simulation
run!(simulation)

# # A simulation of a soda outgassing CO₂ on the counter
# We simulate soda Warming up on the counter using a forcing function
# linear increase from 4°C to 30°C over 12 hours then stops warming
temperature_increase(z, t, p) = ifelse(t <= 12hours, p.∂T∂t * p.Δt , 0.0)

warming = Forcing(
     temperature_increase, 
     parameters=(∂T∂t=1e-6, Δt=10minutes),
     )

# Build the second model for the soda on the counter
model_open_on_the_counter = NonhydrostaticModel(; 
    grid,
    velocities = nothing,
    buoyancy   = nothing,
    closure    = nothing,
    forcing    = (; T=warming),
    tracers    = (:T, :DIC),
    boundary_conditions = (; DIC=dic_bcs),
    )

set!(model_open_on_the_counter, T = Tᵢ, DIC = DICᵢ)

simulation = Simulation(model_open_on_the_counter, Δt=10minutes, stop_time=24hours)

# Add an output writer...
fnm_counter = "soda_on_the_counter.jld2"
outputs  = (; co₂_flux, 
              soda_co₂, 
              atmos_co₂, 
              model_open_on_the_counter.tracers.T, 
              model_open_on_the_counter.tracers.DIC,
              )
simulation.output_writers[:jld2] = JLD2OutputWriter(
    model_open_on_the_counter, outputs;
    filename = fnm_counter,
    schedule = TimeInterval(10minutes),
    overwrite_existing = true)

# ..and don't forget to add a callback to compute the CO₂ flux
add_callback!(simulation, compute_co₂_flux!)

# Run the simulation
run!(simulation)

# ## Visualization
#
# All that's left is to visualize the results.
fridge_soda_co₂ = FieldTimeSeries(fnm_fridge, "soda_co₂")
fridge_atmo_co₂ = FieldTimeSeries(fnm_fridge, "atmos_co₂")
fridge_temp      = FieldTimeSeries(fnm_fridge, "T")

counter_soda_co₂ = FieldTimeSeries(fnm_counter, "soda_co₂")
counter_atmo_co₂ = FieldTimeSeries(fnm_counter, "atmos_co₂")
counter_temp      = FieldTimeSeries(fnm_counter, "T")

t  = fridge_soda_co₂.times
nt = length(t)

fig = Figure(size=(1200, 900))

ax = Axis(fig[1,1], xlabel="Time", ylabel="CO₂ conc [mmol m⁻³]")
lines!(t/(3600), 
       interior(fridge_soda_co₂, 1, 1, 1, :)*1e3; 
       linestyle = :dash, 
       label = "fridge soda",
       )
lines!(t/(3600), 
       interior(counter_soda_co₂, 1, 1, 1, :)*1e3; 
       linestyle = :solid, 
       label = "counter soda",
       )
lines!(t/(3600), 
       interior(fridge_atmo_co₂, 1, 1, 1, :)*1e3;
       linestyle = :dash, 
       label = "fridge saturated")
lines!(t/(3600), 
       interior(counter_atmo_co₂, 1, 1, 1, :)*1e3;
       linestyle = :solid, 
       label = "counter saturated")
axislegend()

ax = Axis(fig[2,1], xlabel="Time", ylabel="Temp (°C)")
lines!(t/(3600), 
       interior(fridge_temp, 1, 1, 1, :),
       linestyle = :dash,
       label = "fridge soda temperature",
       )
lines!(t/(3600), 
       interior(counter_temp, 1, 1, 1, :),
       linestyle = :solid,
       label = "counter soda temperature",
       )
axislegend()
nothing # hide
# The cool soda's CO₂ concentration approaches equilibrium with the atmosphere (the saturated CO₂ concentration) quickly.
# The warming soda continues to outgas, since the solubility of CO₂ decreases with temperature. It'll taste flatter because of the lower CO₂ concentration.
fig
