#using ClimaOceanBiogeochemistry.CarbonSystemSolvers.UniversalRobustCarbonSolver: UniversalRobustCarbonSystem
#using ClimaOceanBiogeochemistry.CarbonSystemSolvers: CarbonSystem, CarbonSystemParameters, CarbonSolverOptions, CarbonCoefficientParameters
using KernelAbstractions: @kernel, @index
using Oceananigans.Utils: launch!
using Oceananigans.Grids: inactive_cell
using Oceananigans.Fields: Field

const exchange_coefficient = 0.337 / 3.6e5 # cm hr⁻¹ / cmhr⁻¹_per_ms⁻¹

"""
    compute_CO₂_flux(
        grid,
        reference_density,
        temperature,
        schmidt_number_dic,
        u_wind_velocity,
        v_wind_velocity,
        surface_wind_speed,
        piston_velocity,
        atmospheric_pCO₂,
        ocean_pCO₂,
        atmospheric_CO₂_solubility,
        oceanic_CO₂_solubility,
        CO₂_flux,
        boundary_condition,
        kˢᶜ = CarbonCoefficientParameters(
            a₀ = 2116.8,
            a₁ = 136.25,
            a₂ = 4.7353,
            a₃ = 9.2307e-2,
            a₄ = 7.555e-4,
            a₅ = 660.0,
            )
    )
(1) Compute the Schmidt number for dissolved inorganic carbon (DIC) based on the temperature Θᶜ (in degrees Celsius).
(2) Compute the wind speed at the ocean surface.
(3) Compute the piston velocity for gas exchange at the ocean surface.
(4) Compute the flux of CO₂ between the atmosphere and the ocean.
(5) Merge the air-sea CO₂ flux into the existing surface DIC flux in-place.

Arguments:
- `grid::RegularCartesianGrid`: The model grid.
- `reference_density::Float64`: The reference density of seawater.
- `temperature::Field{Center, Center, Nothing}`: Temperature in degrees Celsius.
- `kˢᶜ::CarbonCoefficientParameters`: The parameters for the Schmidt number calculation.
- `schmidt_number_dic::Field{Center, Center, Nothing}`: The computed Schmidt number for DIC.
- `u_wind_velocity::Field{Center, Center, Nothing}`: The wind velocity from atmospheric input in the u direction.
- `v_wind_velocity::Field{Center, Center, Nothing}`: The wind velocity from atmospheric input in the v direction.
- `surface_wind_speed::Field{Center, Center, Nothing}`: The wind speed at the ocean surface.
- `piston_velocity::Field{Center, Center, Nothing}`: The piston velocity for gas exchange at the ocean surface.
- `atmospheric_pCO₂::Float64`: The partial pressure of CO₂ in the atmosphere.
- `ocean_pCO₂::Field{Center, Center, Nothing}`: The partial pressure of CO₂ in the ocean.
- `atmospheric_CO₂_solubility::Field{Center, Center, Nothing}`: The solubility of CO₂ in the atmosphere.
- `oceanic_CO₂_solubility::Field{Center, Center, Nothing}`: The solubility of CO₂ in the ocean.
- `CO₂_flux::Field{Center, Center, Nothing}`: The computed flux of CO₂ between the atmosphere and the ocean.
- `boundary_condition::Field{Center, Center, Nothing}`: The boundary condition for the surface DIC flux.
- `kˢᶜ::CarbonCoefficientParameters` (optional): The parameters for the Schmidt number calculation (default values provided).
Notes:
The convention is that a positive flux is upwards (outgassing), and a negative flux is downwards (uptake).
"""
@kernel function compute_CO₂_flux!(
    grid,
    reference_density,
    sea_ice_concentration,
    temperature,
    schmidt_number_dic,
    u_wind_velocity,
    v_wind_velocity,
    surface_wind_speed,
    piston_velocity,
    atmospheric_pCO₂,
    ocean_pCO₂,
    atmospheric_CO₂_solubility,
    oceanic_CO₂_solubility,
    CO₂_flux,
    boundary_condition,
    kˢᶜ = CarbonCoefficientParameters(
            a₀ = 2116.8,
            a₁ = 136.25,
            a₂ = 4.7353,
            a₃ = 9.2307e-2,
            a₄ = 7.555e-4,
            a₅ = 660.0,
        ),
    )

    i, j = @index(Global, NTuple)
    k    = size(grid, 3)
    inactive = inactive_cell(i, j, k, grid)

    # Calculate the Schmidt number for DIC based on the temperature Θᶜ (in degrees Celsius)
    @inbounds begin
        schmidt_number_dic[i, j, 1] = ifelse(
            inactive,
            zero(grid),
            (kˢᶜ.a₀ - 
             kˢᶜ.a₁ * temperature[i, j, k] + 
             kˢᶜ.a₂ * temperature[i, j, k]^2 - 
             kˢᶜ.a₃ * temperature[i, j, k]^3 + 
             kˢᶜ.a₄ * temperature[i, j, k]^4 
            ) /  kˢᶜ.a₅,
        )

        # Compute wind speed at the ocean surface
        surface_wind_speed[i, j, 1] = ifelse(
            inactive,
            zero(grid),
            sqrt(u_wind_velocity[i, j, 1]^2 + 
                 v_wind_velocity[i, j, 1]^2
                ),
        )

        # Compute the piston velocity for gas exchange at the ocean surface
        piston_velocity[i, j, 1] = ifelse(
            inactive,
            zero(grid),
            exchange_coefficient * 
                surface_wind_speed[i, j, 1]^2 / 
                sqrt(schmidt_number_dic[i, j, 1]),
        )

        ## Finally, compute CO₂ flux (-ve for uptake, +ve for outgassing since convention is +ve upwards)
        CO₂_flux[i, j, 1] = ifelse(
            inactive,
            zero(grid),
	        -(one(grid) - sea_ice_concentration[i, j, 1]) * piston_velocity[i, j, 1] * (
                         atmospheric_pCO₂[i,j,1] * atmospheric_CO₂_solubility[i, j, k] - 
                         ocean_pCO₂[i, j, k]     * oceanic_CO₂_solubility[i, j, k]
                      ) * reference_density, # Convert mol kg⁻¹ m s⁻¹ to mol m⁻² s⁻¹
        )

        # Merge air-sea CO₂ flux into the existing surface DIC flux in-place.
        boundary_condition[i, j, 1] += ifelse(
            inactive | isnan(CO₂_flux[i, j, 1]),
            zero(grid),
            CO₂_flux[i, j, 1],
        )
    end
end

"""
    calculate_air_sea_carbon_exchange!(simulation; solver_params = ())
    
Returns the tendency of DIC in the top layer due to air-sea CO₂ flux
using the piston velocity formulation of Wanninkhof (1992) and the
solubility/activity of CO₂ in seawater.
"""
@inline function calculate_air_sea_carbon_exchange!(simulation) #, params = (Pᶜᴼ² = AirSeaCarbonFluxParameters(),)
    grid = simulation.model.ocean.model.grid
    bgc  = simulation.model.ocean.model.biogeochemistry
    ptr  = simulation.model.ocean.model.tracers 
    oce  = simulation.model.interfaces.exchanger.ocean.state
    atm  = simulation.model.interfaces.exchanger.atmosphere.state 
    ice  = simulation.model.interfaces.exchanger.sea_ice.state

    kernel_args = (
        grid,
        bgc.reference_density,
        ice.ℵ,
        oce.T,
        bgc.schmidt_number_dic,
        atm.u,
        atm.v,
        bgc.surface_wind_speed,
        bgc.piston_velocity,
        bgc.atmospheric_pCO₂,
        bgc.ocean_pCO₂,
        bgc.atmospheric_CO₂_solubility,
        bgc.oceanic_CO₂_solubility,
        bgc.CO₂_flux,
        ptr.DIC.boundary_conditions.top.condition,
    )

    launch!(grid.architecture, 
	    grid, 
	    :xy, 
	    compute_CO₂_flux!, 
	    kernel_args...,
    )
    return nothing
end

