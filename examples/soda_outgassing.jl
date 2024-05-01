# # Nutrients, plankton, bacteria, detritus
#
# This example illustrates how to use ClimaOceanBiogeochemistry's
# `CarbonAlkalinityNutrients` model in a single column context.

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

# Build function for CO₂ flux calculation. Dissolved CO₂ in the soda will exchange with the overlying atmosphere
# These are some coefficients and constants that we'll use in the calculation
Base.@kwdef struct Pᶠˡᵘˣᶜᵒ²{FT}
    surface_wind_speed   :: FT = 10. # ms⁻¹
    atmospheric_pCO₂     :: FT = 280e-6 # atm
    exchange_coefficient :: FT = 0.337 # cm hr⁻¹
    salinity             :: FT = 35.0 # psu
    alkalinity           :: FT = 2.35e-3 # mol kg⁻¹
    silicate_conc        :: FT = 0e-6 # mol kg⁻¹
    phosphate_conc       :: FT = 0e-6 # mol kg⁻¹
    initial_pH_guess     :: FT = 8.0
    reference_density    :: FT = 1000. # kg m⁻³
end

"""
    compute_co₂_flux!(simulation)
    
Returns the tendency due to CO₂ flux using the piston velocity 
formulation of Wanninkhof (1992) and the solubility/activity of 
CO₂ that depends on temperature, etc.
"""
@inline function compute_co₂_flux!(simulation)
# conversion factor from cm/hr to m/s
    cmhr⁻¹_per_ms⁻¹ = 1 / 3.6e5

# get coefficients from Pᶠˡᵘˣᶜᵒ² struct
# I really want the option to take these from the model
    (; surface_wind_speed, 
       atmospheric_pCO₂, 
       exchange_coefficient, 
       salinity,
       alkalinity,
       silicate_conc, 
       phosphate_conc,
       initial_pH_guess, 
       reference_density,   
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
    co₂_flux = simulation.model.tracers.DIC.boundary_conditions.top.condition

    # applied pressure in atm (i.e. Pˢᵘʳᶠ-Pᵃᵗᵐ) 
    # Positive when the can is sealed, then zero when the can is opens
    # On average, the 12 ounce soda cans sold in the US tend to have a pressure of roughly 120 kPa when canned at 4 °C
    if iteration(simulation) <= 1
        Δpᵦₐᵣ   = 0.2
    else
        # *pssshhhht* the bottle is opened
        Δpᵦₐᵣ   = 0.0
    end

    # access model fields
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

    # compute sodaic pCO₂ using the UniversalRobustCarbonSystem solver
    # Returns soda pCO₂ (in atm) and atmosphere/soda solubility coefficients (mol kg⁻¹ atm⁻¹)
    (; pCO₂ᵒᶜᵉ, Pᵈⁱᶜₖₛₒₗₐ, Pᵈⁱᶜₖₛₒₗₒ) = 
        UniversalRobustCarbonSystem(
                Θᶜ, Sᴬ, Δpᵦₐᵣ, Cᵀ, Aᵀ, Pᵀ, Siᵀ, pH, pCO₂ᵃᵗᵐ,
                )

    # store the soda and atmospheric pCO₂ into Fields
    soda_pco₂[
        simulation.model.grid.Nx,
        simulation.model.grid.Ny,
        simulation.model.grid.Nz
        ] = pCO₂ᵒᶜᵉ
    atmos_pco₂[
        simulation.model.grid.Nx,
        simulation.model.grid.Ny,
        simulation.model.grid.Nz
        ] = pCO₂ᵃᵗᵐ
    
    # compute schmidt number for DIC
    Sᶜᵈⁱᶜ =  2116.8 - 
              136.25 * Θᶜ + 
                4.7353 * Θᶜ^2 - 
                0.092307 * Θᶜ^3 + 
                7.555e-4 * Θᶜ^4
    
    # compute gas exchange coefficient/piston velocity and correct with Schmidt number
    Kʷ =  (
           (Kʷₐᵥₑ * cmhr⁻¹_per_ms⁻¹) * U₁₀^2
          ) / sqrt(Sᶜᵈⁱᶜ/660)    
    
    # compute co₂ flux (-ve for uptake, +ve for outgassing since convention is +ve upwards in the soda)
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

# ## Model setup represents a 10cm bottle of soda
grid = RectilinearGrid(size = 1,
                       z = (-1meter/10, 0),
                       topology = (Flat, Flat, Bounded))

# For CO₂ exchange, we use a FluxBoundaryCondition for the "top" of the dissolved inorganic carbon (DIC) tracer.
# We'll write a callback later to calculate the flux every time step.
co₂_flux = Field{Center, Center, Nothing}(grid)
dic_bcs  = FieldBoundaryConditions(top = FluxBoundaryCondition(co₂_flux))

### We put the pieces together ###
model_open_in_the_fridge = NonhydrostaticModel(; 
    grid,
    velocities = nothing,
    buoyancy   = nothing,
    closure    = nothing,
    tracers    = (:T, :DIC),
    boundary_conditions = (; DIC=dic_bcs),
    )
# ## Initial conditions
Tᵢ   = 4     # °C
DICᵢ = 2.4   # mol/m³ 

set!(model_open_in_the_fridge, T = Tᵢ, DIC = DICᵢ)

# ## A simulation of a soda outgassing CO₂
simulation = Simulation(model_open_in_the_fridge, Δt=10minutes, stop_time=12hours)

# ## A callback to compute the CO₂ flux
add_callback!(simulation, compute_co₂_flux!)

# These are filled in compute_co₂_flux!
soda_pco₂ = Field{Center, Center, Nothing}(grid)
atmos_pco₂ = Field{Center, Center, Nothing}(grid)

# ## An output writer
fnm_fridge = "soda_in_the_fridge.jld2"
outputs = (; co₂_flux, soda_pco₂, atmos_pco₂, model_open_in_the_fridge.tracers.T)
simulation.output_writers[:jld2] = JLD2OutputWriter(
    model_open_in_the_fridge, outputs;
    filename = fnm_fridge,
    schedule = TimeInterval(10minutes),
    overwrite_existing = true,
    )

run!(simulation)

### We put the pieces together for the second model ###
## Warming of soda temperature
# K s⁻¹ linear increase from 4°C to 30°C over 12 hours
temperature_increase(z, t, p) = p.∂T∂t * p.Δt 

warming = Forcing(
     temperature_increase, 
     parameters=(∂T∂t=1e-6, Δt=10minutes),
     )

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

# ## A simulation of a soda outgassing CO₂
simulation = Simulation(model_open_on_the_counter, Δt=10minutes, stop_time=12hours)

# ## A callback to compute the CO₂ flux
add_callback!(simulation, compute_co₂_flux!)

# These are filled in compute_co₂_flux!
soda_pco₂ = Field{Center, Center, Nothing}(grid)
atmos_pco₂ = Field{Center, Center, Nothing}(grid)

# ## An output writer
fnm_counter = "soda_on_the_counter.jld2"
outputs  = (; co₂_flux, soda_pco₂, atmos_pco₂, model_open_on_the_counter.tracers.T)
simulation.output_writers[:jld2] = JLD2OutputWriter(
    model_open_on_the_counter, outputs;
    filename = fnm_counter,
    schedule = TimeInterval(10minutes),
    overwrite_existing = true)

run!(simulation)

# ## Visualization
#
# All that's left is to visualize the results.
fridge_soda_pco₂ = FieldTimeSeries(fnm_fridge, "soda_pco₂")
fridge_atmo_pco₂ = FieldTimeSeries(fnm_fridge, "atmos_pco₂")
fridge_temp      = FieldTimeSeries(fnm_fridge, "T")

counter_soda_pco₂ = FieldTimeSeries(fnm_counter, "soda_pco₂")
counter_atmo_pco₂ = FieldTimeSeries(fnm_counter, "atmos_pco₂")
counter_temp      = FieldTimeSeries(fnm_counter, "T")

t  = fridge_soda_pco₂.times
nt = length(t)

fig = Figure(size=(1200, 900))

ax = Axis(fig[1,1], xlabel="Time", ylabel="pCO₂ (μatm)")
lines!(t/(3600), 
       interior(fridge_soda_pco₂, 1, 1, 1, :)*1e6; 
       linestyle = :dash, 
       label = "fridge soda pCO₂",
       )
lines!(t/(3600), 
       interior(counter_soda_pco₂, 1, 1, 1, :)*1e6; 
       linestyle = :solid, 
       label = "counter soda pCO₂",
       )
lines!(t/(3600), 
       interior(fridge_atmo_pco₂, 1, 1, 1, :)*1e6; 
       linestyle = :solid, 
       label = "atmos pCO₂")
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

## The cool soda's pCO₂ approaches equilibrium with the atmosphere quickly.
## The the warming soda continues to outgas, since the solubility of CO₂ decreases with temperature.
fig
