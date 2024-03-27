# # Nutrients, plankton, bacteria, detritus
#
# This example illustrates how to use ClimaOceanBiogeochemistry's
# `CarbonAlkalinityNutrients` model in a single column context.

using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients
using ClimaOceanBiogeochemistry.CarbonSystemSolvers.UniversalRobustCarbonSolver: UniversalRobustCarbonSystem
include("../src/carbon_chemistry_coefficients.jl")

using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity

using Printf
using CairoMakie

# ## A single column grid
#
# We set up a single column grid with 4 m grid spacing that's 256 m deep:

grid = RectilinearGrid(size = 64,
                       z = (-256meters, 0),
                       topology = (Flat, Flat, Bounded))


Qᵇ(t) = ifelse(t < 4days, 1e-7, 0.0)
b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵇ))

# ## Buoyancy that depends on temperature and salinity
# We use the `SeawaterBuoyancy` model with a linear equation of state,
# where thermal expansion
αᵀ = 2e-4
#and haline contraction
βˢ = 8e-4

# ## Boundary conditions
#
# We calculate the surface temperature flux associated with surface cooling of
# 200 W m⁻², reference density `ρₒ`, and heat capacity `cᴾ`,
# To illustrate the dynamics of `CarbonAlkalinityNutrients`, this strong convection
# drives turbulent mixing for 4 days, and then abruptly shuts off. 
# Once the convective turbulence dies down, plankton start to grow.
Qʰ = 200.0  # W m⁻², surface _heat_ flux
ρₒ = 1026.0 # kg m⁻³, average density at the surface of the world ocean
cᴾ = 3991.0 # J K⁻¹ kg⁻¹, typical heat capacity for seawater

Qᵀ(t) = ifelse(t < 4days, Qʰ / (ρₒ * cᴾ), 0.0) # K m s⁻¹, surface _temperature_ flux

# Finally, we impose a temperature gradient `dTdz` both initially and at the
# bottom of the domain, culminating in the boundary conditions on temperature,

dTdz = 0.01 # K m⁻¹

T_bcs = FieldBoundaryConditions(top    = FluxBoundaryCondition(Qᵀ),
                                bottom = GradientBoundaryCondition(dTdz))

# For air-sea CO₂ fluxes, we use a FluxBoundaryCondition for the "top" of the DIC tracer.
# We'll write a callback to calculate the flux every time step.
co₂_flux = Field{Center, Center, Nothing}(grid)
top_dic_bc = FluxBoundaryCondition(co₂_flux)
dic_bcs = FieldBoundaryConditions(top=top_dic_bc)

### We put the pieces together ###
# The important line here is `biogeochemistry = CarbonAlkalinityNutrients()`.
# We use all default parameters.
model = HydrostaticFreeSurfaceModel(; grid,
                                    biogeochemistry = CarbonAlkalinityNutrients(; grid),
                                    closure = CATKEVerticalDiffusivity(),
                                    tracers = (:T, :S, :e, :DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe),
                                    buoyancy = SeawaterBuoyancy(
                                        equation_of_state=LinearEquationOfState(
                                            thermal_expansion  = αᵀ,
                                            haline_contraction = βˢ)
                                            ),
                                    boundary_conditions = (; T=T_bcs, DIC=dic_bcs),
                                    )
# ## Initial conditions
#
## Temperature initial condition: a stable density gradient with random noise superposed.
## Random noise damped at top and bottom
Ξ(z) = randn() * z / model.grid.Lz * (1 + z / model.grid.Lz)
Tᵢ(z) = 20 + dTdz * z + dTdz * model.grid.Lz * 1e-6 * Ξ(z)
Sᵢ(z) = 35

# We initialize the model with reasonable nutrient concentrations
DICᵢ(z) = 2.10 # mol/m³ 
ALKᵢ(z) = 2.35 # mol/m³ 
PO₄ᵢ(z) = 2.5e-3 # mol/m³ 
NO₃ᵢ(z) = 36e-3 # mol/m³ 
DOPᵢ(z) = 0.0 # mol/m³ 
POPᵢ(z) = 0.0 # mol/m³ 
Feᵢ(z)  = 0.6e-6 # mol/m³

set!(
    model, 
    T   = Tᵢ,
    S   = Sᵢ, 
    e   = 1e-6, 
    DIC = DICᵢ, 
    ALK = ALKᵢ, 
    PO₄ = PO₄ᵢ, 
    NO₃ = NO₃ᵢ, 
    DOP = DOPᵢ, 
    POP = POPᵢ, 
    Fe  = Feᵢ
    )

# ## A simulation of physical-biological interaction
# 
# We construct a simple simulation that emits a message every 10 iterations
# and outputs tracer fields.

simulation = Simulation(model, Δt=10minutes, stop_time=12days)

progress(sim) = @printf("Iteration: %d, time: %s\n", iteration(sim), prettytime(sim))
simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

# I really want the option to take these from the model
Base.@kwdef struct Pᶠˡᵘˣᶜᵒ²{FT}
    open_water_fraction :: FT = 1.
    surface_wind_speed  :: FT = 5. # ms⁻¹
    atmospheric_pCO₂    :: FT = 280e-6 # atm
    exchange_coefficient:: FT = 0.337 # cm/hr
    silicate_conc       :: FT = 15e-6 # umol/kg
    initial_pH_guess    :: FT = 8.0
end

# Build operation for CO₂ flux calculation, used in callback below to calculate every time step
"""
    compute_co₂_flux!(simulation)
    
Returns the tendency of DIC in the top layer due to air-sea CO₂ flux
using the piston velocity formulation of Wanninkhof (1992) and the
solubility/activity of CO₂ in seawater.
"""
@inline function compute_co₂_flux!(simulation)
# applied pressure in atm (i.e. Pˢᵘʳᶠ-Pᵃᵗᵐ)
    Δpᵦₐᵣ           = 0.0

# conversion factor from cm/hr to m/s
    cmhr⁻¹_per_ms⁻¹ = 1 / 3.6e5

# get coefficients from Pᶠˡᵘˣᶜᵒ² struct
# I really want the option to take these from the model
    (; open_water_fraction, 
       surface_wind_speed, 
       atmospheric_pCO₂, 
       exchange_coefficient, 
       silicate_conc, 
       initial_pH_guess    
       ) = Pᶠˡᵘˣᶜᵒ²()

    fopen   = open_water_fraction 
    U₁₀     = surface_wind_speed  
    pCO₂ᵃᵗᵐ = atmospheric_pCO₂    
    Kʷₐᵥₑ   = exchange_coefficient * cmhr⁻¹_per_ms⁻¹
    Siᵀ     = silicate_conc       
    pH      = initial_pH_guess 

    co₂_flux = model.tracers.DIC.boundary_conditions.top.condition
    flux = co₂_flux.data

    @inbounds for i in 1:co₂_flux.grid.Nx
        @inbounds for j in 1:co₂_flux.grid.Ny
        # access model fields
            Θᶜ       = model.tracers.T[i,j,model.grid.Nz+1]
            Sᴬ       = model.tracers.S[i,j,model.grid.Nz+1]
            Cᵀ       = model.tracers.DIC[i,j,model.grid.Nz+1]
            Aᵀ       = model.tracers.ALK[i,j,model.grid.Nz+1]
            Pᵀ       = model.tracers.PO₄[i,j,model.grid.Nz+1]

        # compute oceanic pCO₂ using the UniversalRobustCarbonSystem solver
            (; pCO₂ᵒᶜᵉ) = 
            UniversalRobustCarbonSystem(
                    Θᶜ, Sᴬ, Δpᵦₐᵣ, Cᵀ, Aᵀ, Pᵀ, Siᵀ, pH, pCO₂ᵃᵗᵐ,
                    )

        # compute schmidt number for DIC
            Sᶜᵈⁱᶜ =  2116.8 - 
                      136.25 * Θᶜ + 
                        4.7353 * Θᶜ^2 - 
                        0.092307 * Θᶜ^3 + 
                        7.555e-4 * Θᶜ^4

        # compute gas exchange coefficient/piston velocity
            Kʷ =  (Kʷₐᵥₑ * U₁₀^2 * fopen ) / sqrt(Sᶜᵈⁱᶜ/660)    
            
            Cᶜᵒᵉᶠᶠ = CarbonChemistryCoefficients(Θᶜ, Sᴬ, Δpᵦₐᵣ)

        # compute co₂ flux
            flux[i,j] = Kʷ * 
                (pCO₂ᵃᵗᵐ * Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖₛₒₗₐ * - 
                 pCO₂ᵒᶜᵉ * Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀ * Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖₛₒₗₒ)
        end
    end

    set!(co₂_flux, flux)
    return nothing
end

add_callback!(simulation, compute_co₂_flux!)

# Output writer
filename = "single_column_carbon_alkalinity_nutrients.jld2"

simulation.output_writers[:fields] = JLD2OutputWriter(model, model.tracers;
                                                      filename,
                                                      schedule = TimeInterval(20minutes),
                                                      overwrite_existing = true)

run!(simulation)

# ## Visualization
#
# All that's left is to visualize the results.

Tt   = FieldTimeSeries(filename, "T")
St   = FieldTimeSeries(filename, "S")
et   = FieldTimeSeries(filename, "e")
DICt = FieldTimeSeries(filename, "DIC")
Alkt = FieldTimeSeries(filename, "ALK")
PO₄t = FieldTimeSeries(filename, "PO₄")
NO₃t = FieldTimeSeries(filename, "NO₃")
DOPt = FieldTimeSeries(filename, "DOP")
POPt = FieldTimeSeries(filename, "POP")
Fet  = FieldTimeSeries(filename, "Fe")

t  = Tt.times
nt = length(t)
z  = znodes(Tt)

fig = Figure(size=(1200, 600))

axb = Axis(fig[1, 1], xlabel="Temperature (K)", ylabel="z (m)")
axe = Axis(fig[1, 2], xlabel="Turbulent kinetic energy (m² s²)")
axP = Axis(fig[1, 3], xlabel="Concentration (mmol)")
axN = Axis(fig[1, 4], xlabel="Nutrient concentration (mmol)")

xlims!(axe, -1e-5, 1e-3)
xlims!(axP, 0, 0.2)

slider = Slider(fig[2, 1:4], range=1:nt, startvalue=1)
n = slider.value

title = @lift @sprintf("Convecting plankton at t = %d hours", t[$n] / hour)
Label(fig[0, 1:4], title)

Tn = @lift interior(Tt[$n], 1, 1, :)
en = @lift interior(et[$n], 1, 1, :)
DICn = @lift interior(DICt[$n], 1, 1, :) 
Alkn = @lift interior(Alkt[$n], 1, 1, :) 
PO₄n = @lift interior(PO₄t[$n], 1, 1, :) 
NO₃n = @lift interior(NO₃t[$n], 1, 1, :) 
DOPn = @lift interior(DOPt[$n], 1, 1, :) 
Fen  = @lift interior(Fet[$n], 1, 1, :)  

lines!(axb, Tn, z)
lines!(axe, en, z)

lines!(axP, DICn, z, label="DIC")
lines!(axP, PO₄n, z, label="Phosphate")
lines!(axP, Fen, z, label="Iron")
axislegend(axP)

lines!(axN, DOPn, z)

record(fig, "carbon_alkalinity_nutrients.mp4", 1:nt, framerate=24) do nn
    n[] = nn
end
nothing #hide

# ![](carbon_alkalinity_nutrients.mp4)
