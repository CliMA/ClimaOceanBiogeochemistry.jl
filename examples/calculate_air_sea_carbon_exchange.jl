using ClimaOceanBiogeochemistry.CarbonSystemSolvers.UniversalRobustCarbonSolver: UniversalRobustCarbonSystem
using ClimaOceanBiogeochemistry.CarbonSystemSolvers: CarbonSystemParameters, CarbonSolverParameters, CarbonCoefficientParameters
using KernelAbstractions: @kernel, @index
using Oceananigans.Utils: launch!
using Oceananigans.Grids: inactive_cell
## Use default carbon coefficient parameters, but pass solver CarbonSolverParameters as an illustration
#   you can pass anything in fieldnames(CarbonSystemParameters), we check it's valid in UniversalRobustCarbonSolver
solver_params = (Sᵒᵖᵗˢ = CarbonSolverParameters(),)

## Supply some coefficients and external data for the CO₂ flux calculation
Base.@kwdef struct AirSeaCarbonFluxParameters
    atmospheric_pCO₂     :: Real
    exchange_coefficient :: Real
    reference_density    :: Real
end
adapt_structure( 
    to, cp::AirSeaCarbonFluxParameters
    ) = AirSeaCarbonFluxParameters(
           adapt(to, cp.atmospheric_pCO₂),
           adapt(to, cp.exchange_coefficient),
           adapt(to, cp.reference_density),
)
@inline function AirSeaCarbonFluxParameters(; 
        atmospheric_pCO₂::Real = 280e-6, # atm
        exchange_coefficient::Real = 0.337 / 3.6e5, # cm hr⁻¹ / cmhr⁻¹_per_ms⁻¹
        reference_density::Real = 1024.5, # kg m⁻³
        )
    return AirSeaCarbonFluxParameters(atmospheric_pCO₂,
                                      exchange_coefficient,
                                      reference_density,
                                      )
end

"""
    compute_schmidt_dic(
        grid,
        schmidt_number_dic, 
        temperature, 
        kˢᶜ
        )

Compute the Schmidt number for dissolved inorganic carbon (DIC) based on the temperature Θᶜ (in degrees Celsius).

Arguments:
- `grid::RectilinearGrid`: The model grid.
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

    while !inactive
        @inbounds schmidt_number_dic[i, j, 1] =  ( 
             kˢᶜ.a₀ - 
             kˢᶜ.a₁ * temperature[i, j, k] + 
             kˢᶜ.a₂ * temperature[i, j, k]^2 - 
             kˢᶜ.a₃ * temperature[i, j, k]^3 + 
             kˢᶜ.a₄ * temperature[i, j, k]^4 
        ) /  kˢᶜ.a₅
    end
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
- `grid::RectilinearGrid`: The model grid.
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
    inactive = inactive_cell(i, j, 1, grid)

    while !inactive
        @inbounds wind_speed[i, j, 1] = sqrt(
            u_wind_velocity[i, j, 1]^2 + 
            v_wind_velocity[i, j, 1]^2
            )
    end
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
    inactive = inactive_cell(i, j, 1, grid)

    while !inactive
        @inbounds piston_velocity[i, j, 1] = exchange_coefficient[i, j, 1] * 
                                   surface_wind_speed[i, j, 1]^2 / 
                                   sqrt(schmidt_number[i, j, 1])
    end
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

Arguments:
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
    Siᵀ,
    pH, 
    atmosphere_pCO₂
    )
    i, j = @index(Global, NTuple)
    k = size(grid, 3)

    inactive = inactive_cell(i, j, k, grid)

    while !inactive
        @inbounds begin
            ## compute oceanic pCO₂ using the UniversalRobustCarbonSystem solver
            CarbonSolved = UniversalRobustCarbonSystem(;
                pH      = pH[i, j, 1], 
                pCO₂ᵃᵗᵐ = atmosphere_pCO₂[i, j, 1],
                Θᶜ      = temperature[i, j, k], 
                Sᴬ      = salinity[i, j, k], 
                Δpᵦₐᵣ   = applied_pressure_bar[i, j, 1], 
                Cᵀ      = DIC[i, j, k]/reference_density, 
                Aᵀ      = ALK[i, j, k]/reference_density, 
                Pᵀ      = PO4[i, j, k]/reference_density, 
                Siᵀ     = Siᵀ[i, j, k]/reference_density,
                solver_params...,
                )

            ocean_pCO₂[i, j, 1]                 = CarbonSolved.pCO₂ᵒᶜᵉ
            atmospheric_CO₂_solubility[i, j, 1] = CarbonSolved.Pᵈⁱᶜₖₛₒₗₐ
            oceanic_CO₂_solubility[i, j, 1]     = CarbonSolved.Pᵈⁱᶜₖₛₒₗₒ
            pH[i, j, 1]                         = CarbonSolved.pH
        end
    end
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

Arguments:
- `grid::RectilinearGrid`: The model grid.
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
    inactive = inactive_cell(i, j, 1, grid)

    while !inactive
    ## compute CO₂ flux (-ve for uptake, +ve for outgassing since convention is +ve upwards)
        @inbounds CO₂_flux[i, j, 1] = - piston_velocity[i, j, 1] * (
                     atmospheric_pCO₂[i, j, 1] * atmospheric_CO₂_solubility[i, j, 1] - 
                     oceanic_pCO₂[i, j, 1]     * oceanic_CO₂_solubility[i, j, 1]
                  ) * reference_density # Convert mol kg⁻¹ m s⁻¹ to mol m⁻² s⁻¹

    end
end

"""
    calculate_air_sea_carbon_exchange!(simulation; solver_params = ())
    
Returns the tendency of DIC in the top layer due to air-sea CO₂ flux
using the piston velocity formulation of Wanninkhof (1992) and the
solubility/activity of CO₂ in seawater.
"""
@inline function calculate_air_sea_carbon_exchange!(simulation, params = (Pᶜᴼ² = AirSeaCarbonFluxParameters(),))
    grid = simulation.model.ocean.model.grid

## get coefficients from CO₂_flux_parameters struct
## I really want the option to take these from the model
    (; atmospheric_pCO₂, 
       exchange_coefficient,
       reference_density,
        ) = params.Pᶜᴼ²

    ## Access model fields
    CO₂_flux = simulation.model.ocean.model.tracers.DIC.boundary_conditions.top.condition
    temperature = simulation.model.ocean.model.tracers.T
    salinity    = simulation.model.ocean.model.tracers.S
    DIC         = simulation.model.ocean.model.tracers.DIC
    ALK         = simulation.model.ocean.model.tracers.ALK
    PO₄         = simulation.model.ocean.model.tracers.PO₄
    
    Siᵀ         = Field{Center, Center, Nothing}(grid)
    Siᵀ         = set!(Siᵀ,
                        [FSiᵀ[i, j, grid.Nz, Time(FSiᵀ.times[1]),
                        ] for i in 1:grid.Nx, 
                              j in 1:grid.Ny],
                )

    applied_pressure_bar = Field{Center, Center, Nothing}(grid)
    kernel_args = (
        grid,
        applied_pressure_bar,
        simulation.model.atmosphere.pressure,
        simulation.model.clock,
        1e-5, # convert Pa to bar
    )
    applied_pressure_bar = applied_pressure_bar - 1
    compute!(applied_pressure_bar)
    
    launch!(arch,
        grid,
        :xy,
        interpolate_atm2oce!,
        kernel_args...,
    )

    u_wind_velocity = Field{Center, Center, Nothing}(grid)     
    kernel_args = (
        grid,
        u_wind_velocity,
        simulation.model.atmosphere.velocities.u,
        Time(simulation.model.clock),
    )

    launch!(arch,
        grid,
        :xy,
        interpolate_atm2oce!,
        kernel_args...,
    )

    v_wind_velocity = Field{Center, Center, Nothing}(grid)     
    kernel_args = (
        grid,
        u_wind_velocity,
        simulation.model.atmosphere.velocities.v,
        Time(simulation.model.clock),
    )

    launch!(arch,
        grid,
        :xy,
        interpolate_atm2oce!,
        kernel_args...,
    )

    ## compute schmidt number for DIC
    schmidt_dic = Field{Center, Center, Nothing}(grid)
    set!(schmidt_dic, 0)

    kernel_args = (
        grid,
        schmidt_dic, 
        temperature,
    )

    launch!(arch, 
	    grid, 
	    :xy, 
	    compute_schmidt_dic!, 
	    kernel_args...,
    )
    
    # Compute wind speed at the ocean surface
    kernel_args = (
        grid,
        u_wind_velocity,
        v_wind_velocity,
        wind_speed,
    )

    launch!(arch,
        grid,
        :xy,
        compute_wind_speed!,
        kernel_args...,
    )

    ## compute gas exchange coefficient/piston velocity and correct with Schmidt number
    average_exchange_coefficient = Field{Center, Center, Nothing}(grid)
    set!(average_exchange_coefficient, exchange_coefficient)

    piston_velocity = Field{Center, Center, Nothing}(grid)
    set!(piston_velocity, 0)

    kernel_args = (
        grid,
        piston_velocity,
        wind_speed,
        average_exchange_coefficient,
        schmidt_dic,
        )

    launch!(arch, 
	    grid,
	    :xy,
	    compute_piston_velocity!, 
	    kernel_args...,
    )

    ## compute oceanic pCO₂ using the UniversalRobustCarbonSystem solver
    set!(atmos_pCO₂, atmospheric_pCO₂)

    #applied_pressure_bar = Field{Center, Center, Nothing}(grid)
    #set!(applied_pressure_bar, applied_pressure)

    atmospheric_CO₂_solubility = Field{Center, Center, Nothing}(grid)
    oceanic_CO₂_solubility     = Field{Center, Center, Nothing}(grid)
    
    # Remove PCO2 from the parameters tuple and pass the rest to the solver
    solver_params = Base.structdiff(params, (Pᶜᴼ²=nothing,))

    kernel_args = (
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
        PO₄, 
        Siᵀ,
        pH, 
        atmos_pCO₂,
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
        atmos_pCO₂, 
        ocean_pCO₂,
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