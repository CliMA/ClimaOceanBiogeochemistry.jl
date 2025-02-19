# # Illustration of using a carbon system solver
#
# This example illustrates how to use ClimaOceanBiogeochemistry's
# `UniversalRobustCarbonSystem` model in a 0-d context.

using ClimaOceanBiogeochemistry.CarbonSystemSolvers.UniversalRobustCarbonSolver: UniversalRobustCarbonSystem
using ClimaOceanBiogeochemistry.CarbonSystemSolvers: CarbonSystemParameters, CarbonSolverParameters, CarbonCoefficientParameters

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
soda_pCO₂ = Field{Center, Center, Nothing}(grid)
atmos_pCO₂ = Field{Center, Center, Nothing}(grid)

# Build function for CO₂ flux calculation. Dissolved CO₂ in the soda will exchange with the overlying atmosphere
# These are some coefficients and constants that we'll use in the calculation
Base.@kwdef struct CO₂_flux_parameters{FT}
    surface_wind_speed   :: FT = 10. # ms⁻¹
    atmospheric_pCO₂     :: FT = 280e-6 # atm
    exchange_coefficient :: FT = 0.337 # cm hr⁻¹
    salinity             :: FT = 0.0 # psu
    alkalinity           :: FT = 2.35e-3 # mol kg⁻¹
    silicate             :: FT = 0e-6 # mol kg⁻¹
    phosphate            :: FT = 0e-6 # mol kg⁻¹
    initial_pH_guess     :: FT = 8.0
    reference_density    :: FT = 1024.5 # kg m⁻³
end

"""
    compute_co₂_flux!(simulation)
    
Returns the tendency due to CO₂ flux using the piston velocity 
formulation of Wanninkhof (1992) and the solubility/activity of 
CO₂ that depends on temperature, etc.
"""
@inline function compute_co₂_flux!(simulation; solver_params = ())
    Nz = size(simulation.model.grid, 3)

## Get coefficients from CO₂_flux_parameters struct
## I really want the option to take these from the model
    (; surface_wind_speed, 
       atmospheric_pCO₂, 
       exchange_coefficient, 
       salinity,
       alkalinity,
       silicate, 
       phosphate,
       initial_pH_guess, 
       reference_density,   
       ) = CO₂_flux_parameters()

    U₁₀      = surface_wind_speed  
    pCO₂ᵃᵗᵐ  = atmospheric_pCO₂    
    Kʷₐᵥₑ    = exchange_coefficient
    Sᴬ       = salinity
    Aᵀ       = alkalinity
    Siᵀ      = silicate   
    Pᵀ       = phosphate  
    pH       = initial_pH_guess
    ρʳᵉᶠ     = reference_density 

    cmhr⁻¹_per_ms⁻¹ = 1 / 3.6e5 # conversion factor from cm/hr to m/s

    co₂_flux = simulation.model.tracers.DIC.boundary_conditions.top.condition
    Θᶜ = simulation.model.tracers.T[1,1,Nz]
    Cᵀ = simulation.model.tracers.DIC[1,1,Nz]/ρʳᵉᶠ

    ## applied pressure in atm (i.e. Pˢᵘʳᶠ-Pᵃᵗᵐ) 
    ## Positive when the can is sealed, then zero when the can is opens
    ## On average, the 12 ounce soda cans sold in the US tend to have a pressure of roughly 120 kPa when canned at 4 °C
    if iteration(simulation) <= 1
        Δpᵦₐᵣ   = 0.2
    else
        ## *pssshhhht* the bottle is opened
        Δpᵦₐᵣ   = 0.0
    end

    ## compute soda pCO₂ using the UniversalRobustCarbonSystem solver
    ## Returns soda pCO₂ (in atm) and atmosphere/soda solubility coefficients (mol kg⁻¹ atm⁻¹)
    (; pCO₂ᵒᶜᵉ, Pᵈⁱᶜₖₛₒₗₐ, Pᵈⁱᶜₖₛₒₗₒ) = UniversalRobustCarbonSystem(
        pH      = pH, 
        pCO₂ᵃᵗᵐ = pCO₂ᵃᵗᵐ,
        Θᶜ      = Θᶜ, 
        Sᴬ      = Sᴬ, 
        Δpᵦₐᵣ   = Δpᵦₐᵣ, 
        Cᵀ      = Cᵀ, 
        Aᵀ      = Aᵀ, 
        Pᵀ      = Pᵀ,
        Siᵀ     = Siᵀ, 
        solver_params...,
    )

    ## store the soda and atmospheric CO₂ concentrations into Fields
    soda_pCO₂[1,1,Nz]  = (pCO₂ᵒᶜᵉ * Pᵈⁱᶜₖₛₒₗₒ ) * ρʳᵉᶠ # Convert mol kg⁻¹ m s⁻¹ to mol m⁻² s⁻¹
    atmos_pCO₂[1,1,Nz] = (pCO₂ᵃᵗᵐ * Pᵈⁱᶜₖₛₒₗₐ) * ρʳᵉᶠ # Convert mol kg⁻¹ to mol m⁻³ 

    ## compute schmidt number for DIC
    kˢᶜ = CarbonCoefficientParameters(
            a₀ = 2116.8,
            a₁ = 136.25,
            a₂ = 4.7353,
            a₃ = 9.2307e-2,
            a₄ = 7.555e-4,
            a₅ = 660.0,
        )

    Cˢᶜᵈⁱᶜ =  ( kˢᶜ.a₀ - 
                kˢᶜ.a₁ * Θᶜ + 
                kˢᶜ.a₂ * Θᶜ^2 - 
                kˢᶜ.a₃ * Θᶜ^3 + 
                kˢᶜ.a₄ * Θᶜ^4 
              )/kˢᶜ.a₅
    
    ## compute gas exchange coefficient/piston velocity and correct with Schmidt number
    Kʷ =  (
           (Kʷₐᵥₑ * cmhr⁻¹_per_ms⁻¹) * U₁₀^2
          ) / sqrt(Cˢᶜᵈⁱᶜ)    
    
    ## compute co₂ flux (-ve for uptake, +ve for outgassing since convention is +ve upwards in the soda)
    co₂_flux[1,1,Nz] = - Kʷ * (
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
             soda_pCO₂, 
             atmos_pCO₂, 
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
              soda_pCO₂, 
              atmos_pCO₂, 
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
fridge_soda_pCO₂ = FieldTimeSeries(fnm_fridge, "soda_pCO₂")
fridge_atmo_co₂ = FieldTimeSeries(fnm_fridge, "atmos_pCO₂")
fridge_temp      = FieldTimeSeries(fnm_fridge, "T")

counter_soda_pCO₂ = FieldTimeSeries(fnm_counter, "soda_pCO₂")
counter_atmo_co₂ = FieldTimeSeries(fnm_counter, "atmos_pCO₂")
counter_temp      = FieldTimeSeries(fnm_counter, "T")

t  = fridge_soda_pCO₂.times
nt = length(t)

fig = Figure(size=(1200, 900))

ax = Axis(fig[1,1], xlabel="Time", ylabel="CO₂ conc [mmol m⁻³]")
lines!(t/(3600), 
       interior(fridge_soda_pCO₂, 1, 1, 1, :)*1e3; 
       linestyle = :dash, 
       label = "fridge soda",
       )
lines!(t/(3600), 
       interior(counter_soda_pCO₂, 1, 1, 1, :)*1e3; 
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
# The cool soda's CO₂ concentration approaches equilibrium with the atmosphere (the saturated CO₂ concentration) quickly.
# The warming soda continues to outgas, since the solubility of CO₂ decreases with temperature. It'll taste flatter because of the lower CO₂ concentration.
#save("soda_outgassing_0d.png", fig)
nothing #hide
fig
# ![](soda_outgassing_0d.png)
