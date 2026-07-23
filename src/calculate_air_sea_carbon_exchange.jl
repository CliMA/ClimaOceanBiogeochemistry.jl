#using ClimaOceanBiogeochemistry.CarbonSystemSolvers.UniversalRobustCarbonSolver: UniversalRobustCarbonSystem
#using ClimaOceanBiogeochemistry.CarbonSystemSolvers: CarbonSystem, CarbonSystemParameters, CarbonSolverOptions, CarbonCoefficientParameters
using KernelAbstractions: @kernel, @index
using Oceananigans.Utils: launch!
using Oceananigans.Grids: inactive_cell
using Oceananigans.Fields: Field

const exchange_coefficient = 0.337 / 3.6e5 # cm hr⁻¹ / cmhr⁻¹_per_ms⁻¹

"""
    compute_schmidt_dic(
        grid,
        schmidt_number_dic, 
        temperature, 
        kˢᶜ
        )

Compute the Schmidt number for dissolved inorganic carbon (DIC) based on the temperature Θᶜ (in degrees Celsius).

Arguments:
- `grid::RegularCartesianGrid`: The model grid.
- `schmidt_number_dic::Field{Center, Center, Nothing}`: The computed Schmidt number for DIC.
- `temperature::Field{Center, Center, Nothing}`: Temperature in degrees Celsius.
- `kˢᶜ::CarbonCoefficientParameters`: The parameters for the Schmidt number calculation.
"""
@kernel function compute_schmidt_dic!(
    grid,
    schmidt_number_dic,
    temperature, 
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
    k = size(grid, 3)
    inactive = inactive_cell(i, j, k, grid)

    @inbounds schmidt_number_dic[i, j, 1] = ifelse(
        inactive,
        zero(grid),
        (kˢᶜ.a₀ - 
         kˢᶜ.a₁ * temperature[i, j, k] + 
         kˢᶜ.a₂ * temperature[i, j, k]^2 - 
         kˢᶜ.a₃ * temperature[i, j, k]^3 + 
         kˢᶜ.a₄ * temperature[i, j, k]^4 
        ) /  kˢᶜ.a₅,
    )
end

"""
    compute_wind_speed(
        grid,
        u_wind_velocity,
        v_wind_velocity,
        wind_speed
        )

Compute the wind speed at the ocean surface.

Arguments:
- `grid::RegularCartesianGrid`: The model grid.
- `u_wind_velocity::Field{Center, Center, Nothing}`: The wind velocity from atmospheric input in the u direction.
- `v_wind_velocity::Field{Center, Center, Nothing}`: The wind velocity from atmospheric input in the v direction.
- `wind_speed::Field{Center, Center, Nothing}`: The wind speed at the ocean surface.
"""
@kernel function compute_wind_speed!(
    grid,
    u_wind_velocity,
    v_wind_velocity,
    wind_speed,
    )

    i, j = @index(Global, NTuple)
    k    = size(grid, 3)
    inactive = inactive_cell(i, j, k, grid)

    @inbounds wind_speed[i, j, 1] = ifelse(
        inactive,
        zero(grid),
        sqrt(u_wind_velocity[i, j, 1]^2 + 
             v_wind_velocity[i, j, 1]^2
            ),
    )
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

Arguments:
- `grid::RegularCartesianGrid`: The model grid.
- `piston_velocity::Field{Center, Center, Nothing}`: The computed piston velocity.
- `surface_wind_speed::Field{Center, Center, Nothing}`: The wind speed at the ocean surface.
- `exchange_coefficient::Float64`: The gas exchange coefficient.
- `schmidt_number::Field{Center, Center, Nothing}`: The Schmidt number.
"""
@kernel function compute_piston_velocity!(
    grid,
    piston_velocity,
    surface_wind_speed,
    schmidt_number,
    )

    i, j = @index(Global, NTuple)
    k    = size(grid, 3)
    inactive = inactive_cell(i, j, k, grid)

    @inbounds piston_velocity[i, j, 1] = ifelse(
        inactive,
        zero(grid),
        exchange_coefficient * 
            surface_wind_speed[i, j, 1]^2 / 
            sqrt(schmidt_number[i, j, 1]),
    )
end

"""
    compute_CO₂_flux(
        CO₂_flux,
        piston_velocity,
        atmospheric_pCO₂,
        oceanic_pCO₂,
        atmospheric_CO₂_solubility,
        oceanic_CO₂_solubility,
        reference_density
    )
        
Compute the flux of CO₂ between the atmosphere and the ocean.

Arguments:
- `reference_density::Float64`: The reference density of seawater.
- `CO₂_flux::Field{Center, Center, Nothing}`: The computed CO₂ flux.
- `piston_velocity::Field{Center, Center, Nothing}`: The piston velocity for gas exchange at the ocean surface.
- `atmospheric_pCO₂::Float64`: The partial pressure of CO₂ in the atmosphere.
- `oceanic_pCO₂::Field{Center, Center, Nothing}`: The partial pressure of CO₂ in the ocean.
- `atmospheric_CO₂_solubility::Field{Center, Center, Nothing}`: The solubility of CO₂ in the atmosphere.
- `oceanic_CO₂_solubility::Field{Center, Center, Nothing}`: The solubility of CO₂ in the ocean.

Notes:
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
    k    = size(grid, 3)
    inactive = inactive_cell(i, j, k, grid)

    ## compute CO₂ flux (-ve for uptake, +ve for outgassing since convention is +ve upwards)
    @inbounds CO₂_flux[i, j, 1] = ifelse(
        inactive,
        zero(grid),
	    -piston_velocity[i, j, 1] * (
                     atmospheric_pCO₂[i,j,1] * atmospheric_CO₂_solubility[i, j, 1] - 
                     oceanic_pCO₂[i, j, 1]     * oceanic_CO₂_solubility[i, j, 1]
                  ) * reference_density, # Convert mol kg⁻¹ m s⁻¹ to mol m⁻² s⁻¹
    )
end

@kernel function combine_dic_fw_and_co2_fluxes!(grid, CO₂_flux, boundary_condition)
    
    i, j = @index(Global, NTuple)
    k    = size(grid, 3)
    inactive = inactive_cell(i, j, k, grid)

    @inbounds boundary_condition[i, j, 1] -= ifelse(
        inactive | isnan(CO₂_flux[i, j, 1]),
        zero(grid),
        CO₂_flux[i, j, 1]
    )
end


"""
    calculate_air_sea_carbon_exchange!(simulation; solver_params = ())
    
Returns the tendency of DIC in the top layer due to air-sea CO₂ flux
using the piston velocity formulation of Wanninkhof (1992) and the
solubility/activity of CO₂ in seawater.
"""
@inline function calculate_air_sea_carbon_exchange!(simulation) #, params = (Pᶜᴼ² = AirSeaCarbonFluxParameters(),)
    grid = simulation.model.ocean.model.grid

### get coefficients from CO₂_flux_parameters struct
### I really want the option to take these from the model
#    (; 
#     # atmospheric_pCO₂, 
#       exchange_coefficient,
#       reference_density,
#        ) = params.Pᶜᴼ²

    ## Access model fields
    #FC          = simulation.model.ocean.model.tracers.DIC.boundary_conditions.top.condition
    #temperature = simulation.model.ocean.model.tracers.T
    #salinity    = simulation.model.ocean.model.tracers.S
    #DIC         = simulation.model.ocean.model.tracers.DIC
    #ALK         = simulation.model.ocean.model.tracers.ALK
    #PO₄         = simulation.model.ocean.model.tracers.PO₄
    
    #Siᵀ         = Field{Center, Center, Center}(grid)
    #update_field_time_series!(FSiᵀ, simulation.model.time)
    #set!(Siᵀ, FSiᵀ)
    
    #pressure_Pa     = simulation.model.fluxes.surface_atmosphere_state.p
    #u_wind_velocity = simulation.model.fluxes.surface_atmosphere_state.u
    #v_wind_velocity = simulation.model.fluxes.surface_atmosphere_state.v
    
    schmidt_dic                = Field{Center, Center, Nothing}(grid)
    piston_velocity            = Field{Center, Center, Nothing}(grid)
    wind_speed                 = Field{Center, Center, Nothing}(grid)

    ## compute schmidt number for DIC
    kernel_args = (
        grid,
        schmidt_dic, 
        simulation.model.ocean.model.tracers.T,
    )

    launch!(grid.architecture, 
	    grid, 
	    :xy, 
	    compute_schmidt_dic!, 
	    kernel_args...,
    )
    
    # Compute wind speed at the ocean surface
    # wind_speed = Field{Center, Center, Nothing}(grid)
    kernel_args = (
        grid,
        simulation.model.interfaces.exchanger.atmosphere.state.u,
        simulation.model.interfaces.exchanger.atmosphere.state.v,
        wind_speed,
    )

    launch!(grid.architecture,
            grid,
            :xy,
            compute_wind_speed!,
            kernel_args...,
    )

    ## compute gas exchange coefficient/piston velocity and correct with Schmidt number
    #piston_velocity = Field{Center, Center, Nothing}(grid)

    kernel_args = (
        grid,
        piston_velocity,
        wind_speed,
        schmidt_dic,
        )

    launch!(grid.architecture, 
	    grid,
	    :xy,
	    compute_piston_velocity!, 
	    kernel_args...,
    )

    ## compute oceanic pCO₂ using the UniversalRobustCarbonSystem solver
    #set!(atmos_pCO₂, atmospheric_pCO₂)

    #applied_pressure_bar = Field{Center, Center, Nothing}(grid)
    #set!(applied_pressure_bar, applied_pressure)

#    # Remove PCO2 from the parameters tuple and pass the rest to the solver
#    #solver_params_only = Base.structdiff(
#				    params, 
#				    (Pᶜᴼ²       = nothing, 
#			             CO₂_flux   = nothing,
#			             ocean_pCO₂ = nothing,
#				     pH         = nothing,),
#				    )
#    kernel_args = (
#        grid,
#        solver_params_only,
#        reference_density,
#        ocean_pCO₂, 
#        atmospheric_CO₂_solubility, 
#        oceanic_CO₂_solubility,
#        simulation.model.ocean.model.tracers.T, 
#        simulation.model.ocean.model.tracers.S, 
#        simulation.model.fluxes.surface_atmosphere_state.p, #applied_pressure_bar, 
#        simulation.model.ocean.model.tracers.DIC, 
#        simulation.model.ocean.model.tracers.ALK,
#        simulation.model.ocean.model.tracers.PO₄, 
#        pH, 
#        atmospheric_pCO₂,
#        )
#
#    launch!(grid.architecture,
#	    grid,
#	    :xy,
#	    solve_ocean_pCO₂!,
#	    kernel_args...,
#    )

    ## compute CO₂ flux (-ve for uptake, +ve for outgassing since convention is +ve upwards)
    kernel_args = (
        grid,
        simulation.model.ocean.model.biogeochemistry.reference_density,
        simulation.model.ocean.model.biogeochemistry.CO₂_flux,
        piston_velocity,
        simulation.model.ocean.model.biogeochemistry.atmospheric_pCO₂, 
        simulation.model.ocean.model.biogeochemistry.ocean_pCO₂,
        simulation.model.ocean.model.biogeochemistry.atmospheric_CO₂_solubility, 
        simulation.model.ocean.model.biogeochemistry.oceanic_CO₂_solubility,
        )

    launch!(grid.architecture, 
	    grid, 
	    :xy, 
	    compute_CO₂_flux!, 
	    kernel_args...,
    )

    # Update surface freshwater forcing on DIC
    kernel_args =(
        grid,
        simulation.model.ocean.model.biogeochemistry.CO₂_flux,
        simulation.model.ocean.model.tracers.DIC.boundary_conditions.top.condition,
    )
    launch!(
        grid.architecture,
        grid,
        :xy,
        combine_dic_fw_and_co2_fluxes!,
        kernel_args...,
    )

    # Done!
    return nothing
end

