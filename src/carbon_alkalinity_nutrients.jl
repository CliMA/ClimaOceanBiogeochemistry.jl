import Oceananigans.Biogeochemistry:
       required_biogeochemical_tracers,
       required_biogeochemical_auxiliary_fields,
       biogeochemical_drift_velocity,
       biogeochemical_auxiliary_fields,
       update_biogeochemical_state!

using Oceananigans.Biogeochemistry: AbstractBiogeochemistry
using Adapt
import Adapt: adapt_structure, adapt
using Oceananigans.BoundaryConditions: ImpenetrableBoundaryCondition, fill_halo_regions!
using Oceananigans.Fields: ConstantField, ZeroField, CenterField
using Oceananigans.Grids: Center, znode, znodes
using Oceananigans.Units: days
using Oceananigans.Architectures: architecture
using ClimaOceanBiogeochemistry.CarbonSystemSolvers.UniversalRobustCarbonSolver: UniversalRobustCarbonSystem
using ClimaOceanBiogeochemistry.CarbonSystemSolvers: CarbonSystem, CarbonSystemParameters, CarbonSolverParameters, CarbonCoefficientParameters
using KernelAbstractions: @kernel, @index
using Oceananigans.Utils: launch!
using Oceananigans.Grids: inactive_cell

#include("calculate_air_sea_carbon_exchange.jl")

const c = Center()
const Pa2bar = 1e-5

struct CarbonAlkalinityNutrients{FT, W, S, M, P} <: AbstractBiogeochemistry
    reference_density                               :: FT # kg m⁻³
    maximum_net_community_production_rate           :: FT # mol PO₄ m⁻³ s⁻¹
    phosphate_half_saturation                       :: FT # mol PO₄ m⁻³
    nitrate_half_saturation                         :: FT # mol NO₃ m⁻³
    iron_half_saturation                            :: FT # mol Fe m⁻³
    incident_PAR                                    :: S # W m⁻²
    PAR_fraction_of_incoming_solar_radiation        :: FT
    PAR_half_saturation                             :: FT  # W m⁻²
    PAR_attenuation_scale                           :: FT  # m
    PAR_percent                                     :: FT
    fraction_of_particulate_export                  :: FT
    dissolved_organic_phosphorus_remin_timescale    :: FT # s⁻¹
    option_of_particulate_remin                     :: FT 
    particulate_organic_phosphorus_remin_timescale  :: FT # s⁻¹
    stoichoimetric_ratio_carbon_to_phosphate        :: FT 
    stoichoimetric_ratio_nitrate_to_phosphate       :: FT 
    stoichoimetric_ratio_phosphate_to_oxygen        :: FT 
    stoichoimetric_ratio_iron_to_phosphate          :: FT 
    stoichoimetric_ratio_carbon_to_nitrate          :: FT 
    stoichoimetric_ratio_carbon_to_oxygen           :: FT 
    stoichoimetric_ratio_carbon_to_iron             :: FT 
    stoichoimetric_ratio_silicate_to_phosphate      :: FT 
    rain_ratio_inorganic_to_organic_carbon          :: FT 
    martin_curve_exponent                           :: FT 
    iron_scavenging_rate                            :: FT # s⁻¹
    ligand_concentration                            :: FT # mol L m⁻³
    ligand_stability_coefficient                    :: FT
    particulate_organic_phosphorus_sinking_velocity :: W  # m s⁻¹
    PAR                                             :: M
    pH                                              :: M
    ocean_pCO₂                                      :: M
    atmospheric_CO₂_solubility                      :: M
    oceanic_CO₂_solubility                          :: M
    carbon_solver_params                            :: P # Named Tuple
    atmospheric_pCO₂                                :: FT # ATM
    CO₂_flux                                        :: M
end

"""
    CarbonAlkalinityNutrients(  reference_density                             = 1024.5,
                                maximum_net_community_production_rate         = 1 / day,
                                phosphate_half_saturation                     = 1e-7 * reference_density,
                                nitrate_half_saturation                       = 1.6e-6 * reference_density,
                                iron_half_saturation                          = 1e-10 * reference_density,
                                incident_PAR                                  = 700.0,
                                PAR_fraction_of_incoming_solar_radiation      = 0.4,
				PAR_half_saturation                           = 10.0,
                                PAR_attenuation_scale                         = 25.0,
                                PAR_percent                                   = 0.01,
                                fraction_of_particulate_export                = 0.33,
                                dissolved_organic_phosphorus_remin_timescale  = 1 / 30day,
                                option_of_particulate_remin                   = 1,
                                particulate_organic_phosphorus_remin_timescale= 0.03 / day, 
                                stoichoimetric_ratio_carbon_to_phosphate      = 106.0
                                stoichoimetric_ratio_nitrate_to_phosphate     = 16.0
                                stoichoimetric_ratio_phosphate_to_oxygen      = 170.0,
                                stoichoimetric_ratio_iron_to_phosphate        = 4.68e-4
                                stoichoimetric_ratio_carbon_to_nitrate        = 106 / 16
                                stoichoimetric_ratio_carbon_to_oxygen         = 106 / 170,
                                stoichoimetric_ratio_carbon_to_iron           = 106 / 1.e-3,
                                stoichoimetric_ratio_silicate_to_phosphate    = 15.0,
                                rain_ratio_inorganic_to_organic_carbon        = 1e-1,
                                martin_curve_exponent                         = 0.84, 
                                iron_scavenging_rate                          = 5e-4 / day,
                                ligand_concentration                          = 1e-9 * reference_density,
                                ligand_stability_coefficient                  = 1e8
                                particulate_organic_phosphorus_sinking_velocity    = -10.0 / day)

Return a seven-tracer biogeochemistry model for the interaction of carbon, alkalinity, and nutrients.

Keyword Arguments
=================

Tracer names
============
* `DIC`: Dissolved Inorganic Carbon

* `ALK`: Alkalinity

* `PO₄`: Phosphate (macronutrient)

* `NO₃`: Nitrate (macronutrient)

* `DOP`: Dissolved Organic Phosphorus

* `POP`: Particulate Organic Phosphorus

* `Fe`: Dissolved Iron (micronutrient)

Biogeochemical functions
========================
* transitions for `DIC`, `ALK`, `PO₄`, `NO₃`, `DOP`, `POP` and `Fe`

* `biogeochemical_drift_velocity` for `POP`, modeling the sinking of detritus at
  a constant `detritus_sinking_speed`.

"""
function CarbonAlkalinityNutrients(; grid,
                                   reference_density                             = 1024.5,
                                   maximum_net_community_production_rate         = 4.e-3 / 365.25days, # mol PO₄ m⁻³ s⁻¹
                                   phosphate_half_saturation                     = 5.e-7 * reference_density, # mol PO₄ m⁻³
                                   nitrate_half_saturation                       = 7.e-6 * reference_density, # mol NO₃ m⁻³
                                   iron_half_saturation                          = 1.e-10 * reference_density, # mol Fe m⁻³
                                   incident_PAR                                  = 700.0, # W m⁻²
                                   PAR_fraction_of_incoming_solar_radiation      = 0.4,
                                   PAR_half_saturation                           = 30.0,  # W m⁻²
                                   PAR_attenuation_scale                         = 25.0,  # m
                                   PAR_percent                                   = 0.01,  # % light: to define the euphotic zone
                                   fraction_of_particulate_export                = 0.33,
                                   dissolved_organic_phosphorus_remin_timescale  = 2. / 365.25days, # s⁻¹
                                   option_of_particulate_remin                   = 1, # POP remin rate: 1 for 1/z, others for constant
                                   particulate_organic_phosphorus_remin_timescale= 0.03 / day, # s⁻¹
                                   stoichoimetric_ratio_carbon_to_phosphate      = 117.0,
                                   stoichoimetric_ratio_nitrate_to_phosphate     = 16.0,
                                   stoichoimetric_ratio_phosphate_to_oxygen      = 170.0, 
                                   stoichoimetric_ratio_iron_to_phosphate        = 4.68e-4,
                                   stoichoimetric_ratio_carbon_to_nitrate        = 117. / 16.,
                                   stoichoimetric_ratio_carbon_to_oxygen         = 117. / 170., 
                                   stoichoimetric_ratio_carbon_to_iron           = 117. / 4.68e-4,
                                   stoichoimetric_ratio_silicate_to_phosphate    = 15.0,
                                   rain_ratio_inorganic_to_organic_carbon        = 1.e-2,
                                   martin_curve_exponent                         = 0.84, 
                                   iron_scavenging_rate                          = 0.2 / 365.25days, # s⁻¹
                                   ligand_concentration                          = 1e-9 * reference_density, # mol L m⁻³
                                   ligand_stability_coefficient                  = 1e8,
                                   particulate_organic_phosphorus_sinking_velocity   = -10.0 / day,
                                   carbon_solver_params                          = (),
                                   atmospheric_pCO₂                              = 280e-6,
				   )
    if incident_PAR isa Number
        surface_PAR = incident_PAR            
        incident_PAR = CenterField(grid)            
        set!(incident_PAR, surface_PAR)            
        fill_halo_regions!(incident_PAR)
    end
    S = typeof(incident_PAR)

    if particulate_organic_phosphorus_sinking_velocity  isa Number
        w₀ = particulate_organic_phosphorus_sinking_velocity 
        no_penetration = ImpenetrableBoundaryCondition()
        bcs = FieldBoundaryConditions(grid, (Center, Center, Face),
                                      top=no_penetration, bottom=no_penetration)
        particulate_organic_phosphorus_sinking_velocity  = ZFaceField(grid, boundary_conditions = bcs)
        set!(particulate_organic_phosphorus_sinking_velocity , w₀)
        fill_halo_regions!(particulate_organic_phosphorus_sinking_velocity )
    end
                            
    FT = eltype(grid)

    pH = CenterField(grid)
    # Set initial first guess of pH
    set!(pH,one(grid)*8)

    ocean_pCO₂                 = CenterField(grid)
    atmospheric_CO₂_solubility = CenterField(grid)
    oceanic_CO₂_solubility     = CenterField(grid)
    PAR                        = CenterField(grid)
    CO₂_flux                   = CenterField(grid)
    M = typeof(ocean_pCO₂)

    return CarbonAlkalinityNutrients(convert(FT, reference_density),
                                     convert(FT, maximum_net_community_production_rate),
                                     convert(FT, phosphate_half_saturation),
                                     convert(FT, nitrate_half_saturation),
                                     convert(FT, iron_half_saturation),
                                     incident_PAR,
                                     convert(FT, PAR_fraction_of_incoming_solar_radiation),
                                     convert(FT, PAR_half_saturation),
                                     convert(FT, PAR_attenuation_scale),     
                                     convert(FT, PAR_percent), 
                                     convert(FT, fraction_of_particulate_export),
                                     convert(FT, dissolved_organic_phosphorus_remin_timescale),
                                     convert(FT, option_of_particulate_remin),
                                     convert(FT, particulate_organic_phosphorus_remin_timescale),
                                     convert(FT, stoichoimetric_ratio_carbon_to_phosphate),
                                     convert(FT, stoichoimetric_ratio_nitrate_to_phosphate),
                                     convert(FT, stoichoimetric_ratio_phosphate_to_oxygen),
                                     convert(FT, stoichoimetric_ratio_iron_to_phosphate),
                                     convert(FT, stoichoimetric_ratio_carbon_to_nitrate),
                                     convert(FT, stoichoimetric_ratio_carbon_to_oxygen),
                                     convert(FT, stoichoimetric_ratio_carbon_to_iron),
                                     convert(FT, stoichoimetric_ratio_silicate_to_phosphate),
                                     convert(FT, rain_ratio_inorganic_to_organic_carbon),
                                     convert(FT, martin_curve_exponent),
                                     convert(FT, iron_scavenging_rate),
                                     convert(FT, ligand_concentration),
                                     convert(FT, ligand_stability_coefficient),
                                     particulate_organic_phosphorus_sinking_velocity,
                                     PAR,
                                     pH,
                                     ocean_pCO₂,
                                     atmospheric_CO₂_solubility,
                                     oceanic_CO₂_solubility,
                                     carbon_solver_params,
                                     convert(FT, atmospheric_pCO₂),
                                     CO₂_flux,
)
end

Adapt.adapt_structure(to, bgc::CarbonAlkalinityNutrients) = 
    CarbonAlkalinityNutrients(adapt(to, bgc.reference_density),
                            adapt(to, bgc.maximum_net_community_production_rate),
                            adapt(to, bgc.phosphate_half_saturation),
                            adapt(to, bgc.nitrate_half_saturation),
                            adapt(to, bgc.iron_half_saturation),
                            adapt(to, bgc.incident_PAR),
                            adapt(to, bgc.PAR_fraction_of_incoming_solar_radiation),
                            adapt(to, bgc.PAR_half_saturation),
                            adapt(to, bgc.PAR_attenuation_scale),
                            adapt(to, bgc.PAR_percent),
                            adapt(to, bgc.fraction_of_particulate_export),
                            adapt(to, bgc.dissolved_organic_phosphorus_remin_timescale),
                            adapt(to, bgc.option_of_particulate_remin),
                            adapt(to, bgc.particulate_organic_phosphorus_remin_timescale),
                            adapt(to, bgc.stoichoimetric_ratio_carbon_to_phosphate),
                            adapt(to, bgc.stoichoimetric_ratio_nitrate_to_phosphate),
                            adapt(to, bgc.stoichoimetric_ratio_phosphate_to_oxygen),
                            adapt(to, bgc.stoichoimetric_ratio_iron_to_phosphate),
                            adapt(to, bgc.stoichoimetric_ratio_carbon_to_nitrate),
                            adapt(to, bgc.stoichoimetric_ratio_carbon_to_oxygen),
                            adapt(to, bgc.stoichoimetric_ratio_carbon_to_iron),
                            adapt(to, bgc.stoichoimetric_ratio_silicate_to_phosphate),
                            adapt(to, bgc.rain_ratio_inorganic_to_organic_carbon),
                            adapt(to, bgc.martin_curve_exponent),
                            adapt(to, bgc.iron_scavenging_rate),
                            adapt(to, bgc.ligand_concentration),
                            adapt(to, bgc.ligand_stability_coefficient),
                            adapt(to, bgc.particulate_organic_phosphorus_sinking_velocity),
                            adapt(to, bgc.PAR),
                            adapt(to, bgc.pH),
                            adapt(to, bgc.ocean_pCO₂),    
                            adapt(to, bgc.atmospheric_CO₂_solubility),
                            adapt(to, bgc.oceanic_CO₂_solubility),
                            adapt(to, bgc.carbon_solver_params),
                            adapt(to, bgc.atmospheric_pCO₂),
                            adapt(to, bgc.CO₂_flux),
)

# Define function methods required by Oceananigans BGC interface
const CAN = CarbonAlkalinityNutrients

"""
Required biogeochemical tracers for the CarbonAlkalinityNutrients model
"""
@inline required_biogeochemical_tracers(::CAN) = (:DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe)

"""
Add a vertical sinking "drift velocity" for the particulate organic phosphate (POP) tracer.
"""
@inline function biogeochemical_drift_velocity(bgc::CAN, ::Val{:POP})
    u = ZeroField()
    v = ZeroField()
    w = bgc.particulate_organic_phosphorus_sinking_velocity 
    return (; u, v, w)
end

"""
Required biogeochemical auxiliary tracers for the CarbonAlkalinityNutrients model
"""
@inline required_biogeochemical_auxiliary_fields(::CAN) = tuple(:PAR, 
								 :pH, 
								 :ocean_pCO₂, 
								 :atmospheric_CO₂_solubility,
								 :oceanic_CO₂_solubility,
								 :CO₂_flux
								)
"""
Set the biogeochemical auxiliary tracers for the CarbonAlkalinityNutrients model
"""
@inline biogeochemical_auxiliary_fields(bgc::CAN) = (; PAR = bgc.PAR,
                                                       pH  = bgc.pH, 
                                                       ocean_pCO₂ = bgc.ocean_pCO₂, 
                                                       atmospheric_CO₂_solubility = bgc.atmospheric_CO₂_solubility,
                                                       oceanic_CO₂_solubility = bgc.oceanic_CO₂_solubility,
                                                       CO₂_flux = bgc.CO₂_flux,
                                                     )

#Functions that fill the biogeochemical auxiliary tracers for the CarbonAlkalinityNutrients model
"""
Update PhotosyntheticallyActiveRatiation field (i.e. light)
"""
@kernel function update_PhotosyntheticallyActiveRatiation!(bgc, PAR⁰, PAR, grid) 
    i, j = @index(Global, NTuple)
    kt   = size(grid, 3)

    λ  = bgc.PAR_attenuation_scale
    zᶜ = znodes(grid, Center(), Center(), Center())
    
    inactive = inactive_cell(i, j, kt, grid)
    @inbounds PAR[i, j, kt] =  ifelse(
		        inactive,
			zero(grid),
    			PAR⁰[i, j, kt] * exp(zᶜ[kt] / λ),
			)

    for k in kt-1:-1:1
        inactive = inactive_cell(i, j, k, grid)
        @inbounds PAR[i, j, k] =  ifelse(
			inactive,
			zero(grid),
			PAR⁰[i, j, kt] * exp(zᶜ[k] / λ),
			)
    end
end 

# Set the PAR to 40% of the downwelling shortwave radiation
@kernel function calculate_par_from_surface_Qs!(grid, parfac, incident_par, Qs)
    i, j = @index(Global, NTuple)
    k    = size(grid, 3)
    inactive = inactive_cell(i, j, k, grid)

    @inbounds incident_par[i, j, k] = ifelse(
        inactive,
        zero(grid),
        parfac * Qs[i, j, 1]
    )
end

# use salt forcing to inform DIC and ALK forcing using surface average values
@kernel function calculate_bgc_fw_forcing!(grid, C₀, A₀, S₀, FS, FC, FA)
    i, j = @index(Global, NTuple)
    k    = size(grid, 3)
    inactive = inactive_cell(i, j, k, grid)

    # Add carbon flux to the CO2 flux
    @inbounds begin
        FC[i, j, 1] = ifelse(
            inactive,
            zero(grid),
            (C₀[i, j, k]/S₀[i, j, k]) * FS[i, j, 1],
        )
        FA[i, j, 1] = ifelse(
            inactive,
            zero(grid),
            (A₀[i, j, k]/S₀[i, j, k]) * FS[i, j, 1],
        )
    end
end

@inline function copy_surface_atmospheric_state_for_bgc!(simulation::Simulation)
    grid = simulation.model.ocean.model.grid

    # Transfer the surface shortwave radiation to the incident PAR field (par is ~40% of Qs)
    kernel_args = (
        grid,
        simulation.model.ocean.model.biogeochemistry.PAR_fraction_of_incoming_solar_radiation,
        simulation.model.ocean.model.biogeochemistry.incident_PAR,
        simulation.model.fluxes.surface_atmosphere_state.Qs,
        )

    launch!(
	    architecture(grid),
        grid,
        :xy,
        calculate_par_from_surface_Qs!,
        kernel_args...
    )

    # Transfer the surface salinity to the surface DIC and ALK fields
    kernel_args =(
        grid,
        simulation.model.ocean.model.tracers.DIC,
        simulation.model.ocean.model.tracers.ALK,
        simulation.model.ocean.model.tracers.S,
        simulation.model.fluxes.total.ocean.tracers.S, # The FW forcing
        simulation.model.ocean.model.tracers.DIC.boundary_conditions.top.condition,
        simulation.model.ocean.model.tracers.ALK.boundary_conditions.top.condition,
    )
    launch!(
        architecture(grid),
        grid,
        :xy,
        calculate_bgc_fw_forcing!,
        kernel_args...,
    )
    return nothing
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
- `grid::RegularCartesianGrid`: The model grid.
- `solver_params::NamedTuple`: The parameters for the UniversalRobustCarbonSystem solver.
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
    pH, 
    atmosphere_pCO₂
    )
    i, j = @index(Global, NTuple)
    k    = size(grid, 3)
    inactive = inactive_cell(i, j, k, grid)

    @inbounds CarbonSolved = ifelse(
        inactive,
        CarbonSystem{eltype(grid)}(
            zero(grid),
            zero(grid),
            zero(grid),
            zero(grid),
            zero(grid),
            zero(grid),
            zero(grid),
            zero(grid),
            zero(grid),
            zero(grid),
            zero(grid),
        ),
	## compute oceanic pCO₂ using the UniversalRobustCarbonSystem solver
        UniversalRobustCarbonSystem(;
            pH      = pH[i, j, k],
            pCO₂ᵃᵗᵐ = atmosphere_pCO₂,
            Θᶜ      = temperature[i, j, k],
            Sᴬ      = salinity[i, j, k],
            Δpᵦₐᵣ   = applied_pressure_bar, #[i, j, k]*Pa2bar,
            Cᵀ      = DIC[i, j, k]/reference_density,
            Aᵀ      = ALK[i, j, k]/reference_density,
            Pᵀ      = PO4[i, j, k]/reference_density,
            Siᵀ     = PO4[i, j, k]*15/reference_density,
            solver_params...,
            ),
    )

    ocean_pCO₂[i, j, k]                 = CarbonSolved.pCO₂ᵒᶜᵉ
    atmospheric_CO₂_solubility[i, j, k] = CarbonSolved.Pᵈⁱᶜₖₛₒₗₐ
    oceanic_CO₂_solubility[i, j, k]     = CarbonSolved.Pᵈⁱᶜₖ₀ #Pᵈⁱᶜₖₛₒₗₒ
    pH[i, j, k]                         = CarbonSolved.pH
end

#@kernel function combine_dic_fw_and_co2_fluxes!(grid, CO₂_flux, boundary_condition)
#    i, j = @index(Global, NTuple)
#    k    = size(grid, 3)
#    inactive = inactive_cell(i, j, k, grid)
#
#    @inbounds boundary_condition[i, j, 1] += ifelse(
#        inactive,
#        zero(grid),
#        CO₂_flux[i, j, k]
#    )
#end

"""
Update BGC auxiliary fields
"""
@inline function update_biogeochemical_state!(bgc::CAN, model)
    arch = architecture(model.grid)
    grid = model.grid
    
    # Calculate depth dependent Photosynthetically Active Radiation
    kernel_args =(
        bgc,
        bgc.incident_PAR,
        bgc.PAR,
        grid,
    )

    launch!(
	arch, 
	grid, 
	:xy, 
	update_PhotosyntheticallyActiveRatiation!,
        kernel_args...,
    )

    # Solve the ocean carbon system and calculate pH, pCO₂, etc
    kernel_args = (
        grid,
        bgc.carbon_solver_params,
        bgc.reference_density,
        bgc.ocean_pCO₂, 
        bgc.atmospheric_CO₂_solubility, 
        bgc.oceanic_CO₂_solubility,
        model.tracers.T, 
        model.tracers.S, 
        zero(grid), #applied_pressure_bar, 
        model.tracers.DIC, 
        model.tracers.ALK,
        model.tracers.PO₄, 
        bgc.pH, 
        bgc.atmospheric_pCO₂,
        )

    launch!(
        arch,
        grid,
        :xy,
        solve_ocean_pCO₂!,
        kernel_args...,
    )

#    # Update surface freshwater forcing on DIC
#    # NOW DONE IN calculate_air_sea_carbon_exchange!
#    kernel_args =(
#        model.grid,
#        model.biogeochemistry.CO₂_flux,
#        model.tracers.DIC.boundary_conditions.top.condition,
#    )
#    launch!(
#        arch,
#        grid,
#        :xy,
#        combine_dic_fw_and_co2_fluxes!,
#        kernel_args...,
#    )
     return nothing
end

# Now define functions that relate to the specific BGC model
#"""
#    PAR(surface_photosynthetically_active_ratiation, 
#        photosynthetically_active_ratiation_attenuation_scale, 
#        depth)
#
#Calculate the photosynthetically active radiation (PAR) at a given depth due to attenuation.
#"""
#@inline function PAR(surface_photosynthetically_active_ratiation, 
#                     photosynthetically_active_ratiation_attenuation_scale, 
#                     depth)
#
#    I₀ = surface_photosynthetically_active_ratiation
#    λ  = photosynthetically_active_ratiation_attenuation_scale
#    z  = depth
#    return I₀ * exp(z / λ)
#end

"""
    net_community_production(grid,
			     maximum_net_community_production_rate,
                             light_half_saturation, 
                             phosphate_half_saturation, 
                             nitrate_half_saturation, 
                             iron_half_saturation, 
                             photosynthetically_active_radiation, 
                             phosphate_concentration, 
                             nitrate_concentration, 
                             iron_concentration)

"""
@inline function net_community_production(grid,
					  maximum_net_community_production_rate,
                                          light_half_saturation, 
                                          phosphate_half_saturation, 
                                          nitrate_half_saturation, 
                                          iron_half_saturation, 
                                          photosynthetically_active_radiation, 
                                          phosphate_concentration, 
                                          nitrate_concentration, 
                                          iron_concentration)
    μᵖ=maximum_net_community_production_rate
    kᴵ=light_half_saturation
    kᴾ=phosphate_half_saturation
    kᴺ=nitrate_half_saturation
    kᶠ=iron_half_saturation
    I =photosynthetically_active_radiation
    PO₄ =phosphate_concentration
    NO₃ =nitrate_concentration
    Feₜ =iron_concentration

    # First, adjust the concentrations to be non-zero
    I_nonzero = max(zero(grid), I)
    P_nonzero = max(zero(grid), PO₄)
    N_nonzero = max(zero(grid), NO₃)
    F_nonzero = max(zero(grid), Feₜ)

    # calculate the limitation terms
    Lₗᵢₘ = I_nonzero / (I_nonzero + kᴵ)
    Pₗᵢₘ = P_nonzero / (P_nonzero + kᴾ)
    Nₗᵢₘ = N_nonzero / (N_nonzero + kᴺ)
    Fₗᵢₘ = F_nonzero / (F_nonzero + kᶠ)

    # return the net community production
    return max(zero(grid),
        μᵖ * Lₗᵢₘ * min(Pₗᵢₘ, Nₗᵢₘ, Fₗᵢₘ)
    )
end

"""
    dissolved_organic_phosphorus_remin(remineralization_rate, 
                                      dissolved_organic_phosphorus_concentration)
Calculate the remineralization of dissolved organic phosphorus.
"""
@inline dissolved_organic_phosphorus_remin(remineralization_rate, 
                                          dissolved_organic_phosphorus_concentration) = 
        max(0, remineralization_rate * dissolved_organic_phosphorus_concentration)

"""
Calculate remineralization of particulate organic phosphorus according to 
    1) rate decreases with depth (1/z+z₀), 
or  2) a first-order rate constant.
"""
@inline function particulate_organic_phosphorus_remin(grid,
						      option_of_particulate_remin,
                                                      particulate_organic_phosphorus_remin_timescale,  
                                                      martin_curve_exponent,
                                                      particulate_organic_phosphorus_sinking_velocity,
                                                      PAR_attenuation_scale,
                                                      depth, percent_light,
                                                      particulate_organic_phosphorus_concentration)

        Rᵣ = option_of_particulate_remin                                  
        r = particulate_organic_phosphorus_remin_timescale
        b = martin_curve_exponent    
        wₛ = particulate_organic_phosphorus_sinking_velocity
        λ = PAR_attenuation_scale
        z  = depth
        fᵢ= percent_light
        z₀ = log(fᵢ)*λ # The base of the euphotic layer depth (z₀) where PAR is degraded down to 1%     
        POP = particulate_organic_phosphorus_concentration

    return ifelse(Rᵣ == one(grid), max(zero(grid), b * wₛ / (z + z₀) * POP), max(zero(grid), r * POP))
end


"""
Calculate remineralization of particulate inorganic carbon.
"""
@inline particulate_inorganic_carbon_remin() = 0.0

"""
    iron_scavenging(iron_scavenging_rate, 
                    iron_concentration, 
                    ligand_concentration, 
                    ligand_stability_coefficient)

Calculate the scavenging loss of iron. Iron scavenging depends on free iron, which 
involves solving a quadratic equation in terms of ligand concentration and stability 
coefficient. Ligand-complexed iron is protected from being scavenged.
"""
@inline function iron_scavenging(iron_scavenging_rate, 
                                 iron_concentration, 
                                 ligand_concentration, 
                                 ligand_stability_coefficient) 
    kˢᶜᵃᵛ = iron_scavenging_rate
    Fe    = iron_concentration
    Lᶠᵉ   = ligand_concentration
    β     = ligand_stability_coefficient

    # solve for the equilibrium free iron concentration
       # β = FeL / (Feᶠʳᵉᵉ * Lᶠʳᵉᵉ)
       # Lᶠᵉ = FeL + Lᶠʳᵉᵉ
       # Fe = FeL + Feᶠʳᵉᵉ
       # --> R₁(Feᶠʳᵉᵉ)² + R₂ Feᶠʳᵉᵉ + R₃ = 0
       β⁻¹ = 1/β
       R₁  = 1
       R₂  =  (Lᶠᵉ + β⁻¹ - Fe) 
       R₃  = -(Fe * β⁻¹) 

       # simple quadratic solution for roots
       discriminant = sqrt( R₂*R₂ - ( 4*R₁*R₃ ))

       # directly solve for the free iron concentration
       Feᶠʳᵉᵉ = (-R₂ + discriminant) / (2*R₁) 

       # return the linear scavenging rate (net scavenging)
       return (kˢᶜᵃᵛ * Feᶠʳᵉᵉ)
end

"""
Add surface input of iron. This sould be a boundary condition, but for now we just add a constant source.
"""
@inline iron_sources() = 0.0

"""
Tracer sources and sinks for dissolved inorganic carbon (DIC)
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:DIC}, clock, fields)
    μᵖ = bgc.maximum_net_community_production_rate
    kᴺ = bgc.phosphate_half_saturation
    kᴾ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    kᴵ = bgc.PAR_half_saturation
    I₀ = bgc.incident_PAR
    λ = bgc.PAR_attenuation_scale
    fᵢ= bgc.PAR_percent
    γ = bgc.dissolved_organic_phosphorus_remin_timescale
    r = bgc.particulate_organic_phosphorus_remin_timescale
    α = bgc.fraction_of_particulate_export
    Rᶜᴾ = bgc.stoichoimetric_ratio_carbon_to_phosphate   
    Rᴺᴾ = bgc.stoichoimetric_ratio_nitrate_to_phosphate  
    Rᴾᴼ = bgc.stoichoimetric_ratio_phosphate_to_oxygen   
    Rᶜᴺ = bgc.stoichoimetric_ratio_carbon_to_nitrate     
    Rᶜᴼ = bgc.stoichoimetric_ratio_carbon_to_oxygen      
    Rᶜᶠ = bgc.stoichoimetric_ratio_carbon_to_iron        
    Rᶜᵃᶜᵒ³ = bgc.rain_ratio_inorganic_to_organic_carbon 
    b = bgc.martin_curve_exponent    
    wₛ = bgc.particulate_organic_phosphorus_sinking_velocity
    Rᵣ = bgc.option_of_particulate_remin    

    @inbounds begin
       z   = znode(i, j, k, grid, c, c, c)
       # Available photosynthetic radiation
       I   = bgc.PAR[i, j, k]
       PO₄ = fields.PO₄[i, j, k]
       NO₃ = fields.NO₃[i, j, k]
       Feₜ = fields.Fe[i, j, k]
       DOP = fields.DOP[i, j, k]
       POP = fields.POP[i, j, k]
       Wₛ  = wₛ[i,j,k]
    end

    return Rᶜᴾ * (dissolved_organic_phosphorus_remin(γ, DOP) -
                 (1 + Rᶜᵃᶜᵒ³ * α) * net_community_production(grid, μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, PO₄, NO₃, Feₜ) +
                 particulate_organic_phosphorus_remin(grid, Rᵣ, r, b, Wₛ, λ, z, fᵢ, POP)) +
           particulate_inorganic_carbon_remin()
end

"""
Tracer sources and sinks for alkalinity (ALK)
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:ALK}, clock, fields)
    μᵖ = bgc.maximum_net_community_production_rate
    kᴺ = bgc.phosphate_half_saturation
    kᴾ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    kᴵ = bgc.PAR_half_saturation
    I₀ = bgc.incident_PAR
    λ = bgc.PAR_attenuation_scale
    fᵢ= bgc.PAR_percent
    γ = bgc.dissolved_organic_phosphorus_remin_timescale
    r = bgc.particulate_organic_phosphorus_remin_timescale
    α = bgc.fraction_of_particulate_export
    Rᶜᴾ = bgc.stoichoimetric_ratio_carbon_to_phosphate   
    Rᴺᴾ = bgc.stoichoimetric_ratio_nitrate_to_phosphate  
    Rᴾᴼ = bgc.stoichoimetric_ratio_phosphate_to_oxygen   
    Rᶜᴺ = bgc.stoichoimetric_ratio_carbon_to_nitrate     
    Rᶜᴼ = bgc.stoichoimetric_ratio_carbon_to_oxygen      
    Rᶜᶠ = bgc.stoichoimetric_ratio_carbon_to_iron        
    Rᶜᵃᶜᵒ³ = bgc.rain_ratio_inorganic_to_organic_carbon     
    b = bgc.martin_curve_exponent    
    wₛ = bgc.particulate_organic_phosphorus_sinking_velocity
    Rᵣ = bgc.option_of_particulate_remin

    @inbounds begin
       z = znode(i, j, k, grid, c, c, c)
       # Available photosynthetic radiation
       I   = bgc.PAR[i, j, k]
       PO₄ = fields.PO₄[i, j, k]
       NO₃ = fields.NO₃[i, j, k]
       Feₜ = fields.Fe[i, j, k]
       DOP = fields.DOP[i, j, k]
       POP = fields.POP[i, j, k]
       Wₛ  = wₛ[i,j,k]
    end

    return -Rᴺᴾ * (
        - (1 + Rᶜᵃᶜᵒ³ * α) * net_community_production(grid, μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, PO₄, NO₃, Feₜ) +
        dissolved_organic_phosphorus_remin(γ, DOP) +
        particulate_organic_phosphorus_remin(grid, Rᵣ, r, b, Wₛ, λ, z, fᵢ, POP)
        ) + 2 * particulate_inorganic_carbon_remin()
end

"""
Tracer sources and sinks for dissolved inorganic phosphate (PO₄)
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:PO₄}, clock, fields)
    μᵖ = bgc.maximum_net_community_production_rate
    kᴺ = bgc.phosphate_half_saturation
    kᴾ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    kᴵ = bgc.PAR_half_saturation
    I₀ = bgc.incident_PAR
    λ = bgc.PAR_attenuation_scale
    fᵢ= bgc.PAR_percent
    γ = bgc.dissolved_organic_phosphorus_remin_timescale
    r = bgc.particulate_organic_phosphorus_remin_timescale
    α = bgc.fraction_of_particulate_export
    Rᶜᴾ = bgc.stoichoimetric_ratio_carbon_to_phosphate   
    Rᴺᴾ = bgc.stoichoimetric_ratio_nitrate_to_phosphate  
    Rᴾᴼ = bgc.stoichoimetric_ratio_phosphate_to_oxygen   
    Rᶜᴺ = bgc.stoichoimetric_ratio_carbon_to_nitrate     
    Rᶜᴼ = bgc.stoichoimetric_ratio_carbon_to_oxygen      
    Rᶜᶠ = bgc.stoichoimetric_ratio_carbon_to_iron        
    Rᶜᵃᶜᵒ³ = bgc.rain_ratio_inorganic_to_organic_carbon 
    b = bgc.martin_curve_exponent 
    wₛ = bgc.particulate_organic_phosphorus_sinking_velocity   
    Rᵣ = bgc.option_of_particulate_remin    

    @inbounds begin
       z = znode(i, j, k, grid, c, c, c)
       # Available photosynthetic radiation
       I   = bgc.PAR[i, j, k]
       PO₄ = fields.PO₄[i, j, k]
       NO₃ = fields.NO₃[i, j, k]
       Feₜ = fields.Fe[i, j, k]
       DOP = fields.DOP[i, j, k]
       POP = fields.POP[i, j, k]
       Wₛ  = wₛ[i,j,k]
    end

    return - net_community_production(grid, μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, PO₄, NO₃, Feₜ) +
           dissolved_organic_phosphorus_remin(γ, DOP) +
           particulate_organic_phosphorus_remin(grid, Rᵣ, r, b, Wₛ, λ, z, fᵢ, POP)
end

"""
Tracer sources and sinks for dissolved inorganic nitrate (NO₃)
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:NO₃}, clock, fields)
    μᵖ = bgc.maximum_net_community_production_rate
    kᴺ = bgc.phosphate_half_saturation
    kᴾ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    kᴵ = bgc.PAR_half_saturation
    I₀ = bgc.incident_PAR
    λ = bgc.PAR_attenuation_scale
    fᵢ= bgc.PAR_percent
    γ = bgc.dissolved_organic_phosphorus_remin_timescale
    r = bgc.particulate_organic_phosphorus_remin_timescale
    α = bgc.fraction_of_particulate_export
    Rᶜᴾ = bgc.stoichoimetric_ratio_carbon_to_phosphate   
    Rᴺᴾ = bgc.stoichoimetric_ratio_nitrate_to_phosphate  
    Rᴾᴼ = bgc.stoichoimetric_ratio_phosphate_to_oxygen   
    Rᶜᴺ = bgc.stoichoimetric_ratio_carbon_to_nitrate     
    Rᶜᴼ = bgc.stoichoimetric_ratio_carbon_to_oxygen      
    Rᶜᶠ = bgc.stoichoimetric_ratio_carbon_to_iron        
    Rᶜᵃᶜᵒ³ = bgc.rain_ratio_inorganic_to_organic_carbon  
    b = bgc.martin_curve_exponent    
    wₛ = bgc.particulate_organic_phosphorus_sinking_velocity
    Rᵣ = bgc.option_of_particulate_remin   

    @inbounds begin
       z = znode(i, j, k, grid, c, c, c)
       # Available photosynthetic radiation
       I   = bgc.PAR[i, j, k]
       PO₄ = fields.PO₄[i, j, k]
       NO₃ = fields.NO₃[i, j, k]
       Feₜ = fields.Fe[i, j, k]
       DOP = fields.DOP[i, j, k]
       POP = fields.POP[i, j, k]
       Wₛ  = wₛ[i,j,k]
    end

    return Rᴺᴾ * (
           - net_community_production(grid, μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, PO₄, NO₃, Feₜ) +
           dissolved_organic_phosphorus_remin(γ, DOP) +
           particulate_organic_phosphorus_remin(grid, Rᵣ, r, b, Wₛ, λ, z, fᵢ, POP))
end

"""
Tracer sources and sinks for dissolved iron (FeT)
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:Fe}, clock, fields)
    μᵖ = bgc.maximum_net_community_production_rate
    kᴺ = bgc.phosphate_half_saturation
    kᴾ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    kᴵ = bgc.PAR_half_saturation
    I₀ = bgc.incident_PAR
    λ  = bgc.PAR_attenuation_scale
    fᵢ= bgc.PAR_percent
    γ  = bgc.dissolved_organic_phosphorus_remin_timescale
    α  = bgc.fraction_of_particulate_export
    γ = bgc.dissolved_organic_phosphorus_remin_timescale
    r = bgc.particulate_organic_phosphorus_remin_timescale
    Rᶜᴾ = bgc.stoichoimetric_ratio_carbon_to_phosphate   
    Rᴺᴾ = bgc.stoichoimetric_ratio_nitrate_to_phosphate  
    Rᴾᴼ = bgc.stoichoimetric_ratio_phosphate_to_oxygen   
    Rᶜᴺ = bgc.stoichoimetric_ratio_carbon_to_nitrate     
    Rᶜᴼ = bgc.stoichoimetric_ratio_carbon_to_oxygen      
    Rᶜᶠ = bgc.stoichoimetric_ratio_carbon_to_iron
    Rᶠᴾ = bgc.stoichoimetric_ratio_iron_to_phosphate
    Lᶠᵉ   = bgc.ligand_concentration
    β     = bgc.ligand_stability_coefficient
    kˢᶜᵃᵛ = bgc.iron_scavenging_rate
    b = bgc.martin_curve_exponent    
    wₛ = bgc.particulate_organic_phosphorus_sinking_velocity
    Rᵣ = bgc.option_of_particulate_remin

    @inbounds begin
       z = znode(i, j, k, grid, c, c, c)
       # Available photosynthetic radiation
       I   = bgc.PAR[i, j, k]
       PO₄ = fields.PO₄[i, j, k]
       NO₃ = fields.NO₃[i, j, k]
       Feₜ = fields.Fe[i, j, k]
       DOP = fields.DOP[i, j, k]
       POP = fields.POP[i, j, k]
       Wₛ  = wₛ[i,j,k]
    end

    return Rᶠᴾ * (
                - net_community_production(grid, μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, PO₄, NO₃, Feₜ) +
                  dissolved_organic_phosphorus_remin(γ, DOP) + 
                  particulate_organic_phosphorus_remin(grid, Rᵣ, r, b, Wₛ, λ, z, fᵢ, POP)
                 ) + 
                 iron_sources() -
                 iron_scavenging(kˢᶜᵃᵛ, Feₜ, Lᶠᵉ, β)
end


"""
Tracer sources and sinks for dissolved organic phosphorus (DOP)
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:DOP}, clock, fields)
    μᵖ = bgc.maximum_net_community_production_rate
    kᴺ = bgc.phosphate_half_saturation
    kᴾ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    kᴵ = bgc.PAR_half_saturation
    I₀ = bgc.incident_PAR
    λ = bgc.PAR_attenuation_scale
    γ = bgc.dissolved_organic_phosphorus_remin_timescale
    α = bgc.fraction_of_particulate_export

   @inbounds begin
       z = znode(i, j, k, grid, c, c, c)
       # Available photosynthetic radiation
       I   = bgc.PAR[i, j, k]
       PO₄ = fields.PO₄[i, j, k]
       NO₃ = fields.NO₃[i, j, k]
       Feₜ = fields.Fe[i, j, k]
       DOP = fields.DOP[i, j, k]
       POP = fields.POP[i, j, k]
    end    

    return - dissolved_organic_phosphorus_remin(γ, DOP) +
             (1 - α) * net_community_production(grid, μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, PO₄, NO₃, Feₜ)
end

"""
Tracer sources and sinks for Particulate Organic Phosphorus (POP).
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:POP}, clock, fields)
    μᵖ = bgc.maximum_net_community_production_rate
    kᴺ = bgc.phosphate_half_saturation
    kᴾ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    kᴵ = bgc.PAR_half_saturation
    I₀ = bgc.incident_PAR
    λ = bgc.PAR_attenuation_scale
    fᵢ= bgc.PAR_percent
    r = bgc.particulate_organic_phosphorus_remin_timescale
    α = bgc.fraction_of_particulate_export   
    b = bgc.martin_curve_exponent      
    wₛ = bgc.particulate_organic_phosphorus_sinking_velocity
    Rᵣ = bgc.option_of_particulate_remin
    
    @inbounds begin
       z = znode(i, j, k, grid, c, c, c)
       # Available photosynthetic radiation
       I   = bgc.PAR[i, j, k]
       PO₄ = fields.PO₄[i, j, k]
       NO₃ = fields.NO₃[i, j, k]
       Feₜ = fields.Fe[i, j, k]
       DOP = fields.DOP[i, j, k]
       POP = fields.POP[i, j, k]
       Wₛ  = wₛ[i,j,k]
    end
    return α * net_community_production(grid, μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, PO₄, NO₃, Feₜ) -
           particulate_organic_phosphorus_remin(grid, Rᵣ, r, b, Wₛ, λ, z, fᵢ, POP)
end



