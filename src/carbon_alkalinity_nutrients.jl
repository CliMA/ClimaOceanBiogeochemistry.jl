import Oceananigans.Biogeochemistry:
       required_biogeochemical_tracers,
       required_biogeochemical_auxiliary_fields,
       biogeochemical_drift_velocity,
       biogeochemical_auxiliary_fields,
       update_biogeochemical_state!,
       update_tendencies!

using Oceananigans.Biogeochemistry: AbstractBiogeochemistry
using Adapt
import Adapt: adapt_structure, adapt
using Oceananigans.BoundaryConditions: fill_halo_regions!, ImpenetrableBoundaryCondition
using Oceananigans.Fields: ConstantField, ZeroField, CenterField
using Oceananigans.Grids: Center, znode, znodes
using Oceananigans.Units: days
using Oceananigans.Architectures: architecture
using ClimaOceanBiogeochemistry.CarbonSystemSolvers.UniversalRobustCarbonSolver: UniversalRobustCarbonSystem
using ClimaOceanBiogeochemistry.CarbonSystemSolvers: CarbonSystem, CarbonSystemParameters, CarbonSolverOptions, CarbonCoefficientParameters, CarbonChemistryCoefficients
using KernelAbstractions: @kernel, @index
using Oceananigans.Utils: launch!
using Oceananigans.Grids: inactive_cell

#include("calculate_air_sea_carbon_exchange.jl")

const Pa2bar = 1e-5

struct CarbonAlkalinityNutrients{FT, W, S, D, P} <: AbstractBiogeochemistry
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
    particulate_inorganic_carbon_sinking_velocity   :: W  # m s⁻¹
    calcite_dissolution_rate_saturation_threshold   :: FT
    calcite_dissolution_rate_above_saturation       :: FT
    calcite_dissolution_rate_below_saturation       :: FT
    calcite_dissolution_exponent_above_saturation   :: FT
    calcite_dissolution_exponent_below_saturation   :: FT 
    calcite_saturation_state                        :: D
    PAR                                             :: D
    NCP                                             :: D
    ocean_pH                                        :: D
    ocean_pCO₂                                      :: D
    atmospheric_CO₂_solubility                      :: D
    oceanic_CO₂_solubility                          :: D
    carbon_system_params                            :: P # Named Tuple
    schmidt_number_dic                              :: S
    piston_velocity                                 :: S
    surface_wind_speed                              :: S
    atmospheric_pCO₂                                :: S # ATM
    CO₂_flux                                        :: S
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

Return an eight-tracer biogeochemistry model for the interaction of carbon,
alkalinity, nutrients, and particulate inorganic carbon (PIC).

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

* `PIC`: Particulate Inorganic Carbon

* `Fe`: Dissolved Iron (micronutrient)

Biogeochemical functions
========================
* transitions for `DIC`, `ALK`, `PO₄`, `NO₃`, `DOP`, `POP`, `PIC`, and `Fe`

* `biogeochemical_drift_velocity` for `POP`, modeling the sinking of detritus at
  a constant `detritus_sinking_speed`. Same for PIC

"""
function CarbonAlkalinityNutrients(; grid,
                                   reference_density                             = 1024.5,
                                   maximum_net_community_production_rate         = 4.e-3 / 365.25days, # mol PO₄ m⁻³ s⁻¹
                                   phosphate_half_saturation                     = 5.e-7 * reference_density, # mol PO₄ m⁻³
                                   nitrate_half_saturation                       = 8.e-6 * reference_density, # mol NO₃ m⁻³
                                   iron_half_saturation                          = 1.e-10 * reference_density, # mol Fe m⁻³
                                   incident_PAR                                  = 700.0, # W m⁻²
                                   PAR_fraction_of_incoming_solar_radiation      = 0.4,
                                   PAR_half_saturation                           = 30.0,  # W m⁻²
                                   PAR_attenuation_scale                         = 25.0,  # m
                                   PAR_percent                                   = 0.01,  # % light: to define the euphotic zone
                                   fraction_of_particulate_export                = 0.33,
                                   dissolved_organic_phosphorus_remin_timescale  = 2. / 365.25days, # s⁻¹
                                   option_of_particulate_remin                   = 1, # POP remin rate: 1 for 1/z, others for constant
                                   particulate_organic_phosphorus_remin_timescale= 0.03 / days, # s⁻¹
                                   stoichoimetric_ratio_carbon_to_phosphate      = 117.0,
                                   stoichoimetric_ratio_nitrate_to_phosphate     = 16.0,
                                   stoichoimetric_ratio_phosphate_to_oxygen      = 170.0, 
                                   stoichoimetric_ratio_iron_to_phosphate        = 4.68e-4,
                                   stoichoimetric_ratio_carbon_to_nitrate        = 117. / 16.,
                                   stoichoimetric_ratio_carbon_to_oxygen         = 117. / 170., 
                                   stoichoimetric_ratio_carbon_to_iron           = 117. / 4.68e-4,
                                   stoichoimetric_ratio_silicate_to_phosphate    = 15.0,
                                   rain_ratio_inorganic_to_organic_carbon        = 1.e-1,
                                   martin_curve_exponent                         = 0.84, 
                                   iron_scavenging_rate                          = 0.2 / 365.25days, # s⁻¹
                                   ligand_concentration                          = 1e-9 * reference_density, # mol L m⁻³
                                   ligand_stability_coefficient                  = 1e8,
                                   particulate_organic_phosphorus_sinking_velocity = -10.0 / days,
                                   particulate_inorganic_carbon_sinking_velocity = -10.0 / days,
                                   calcite_dissolution_rate_saturation_threshold = 0.8272,
                                   calcite_dissolution_rate_above_saturation     = 5.22e-9,
                                   calcite_dissolution_rate_below_saturation     = 1.65e-5,
                                   calcite_dissolution_exponent_above_saturation = 0.11,
                                   calcite_dissolution_exponent_below_saturation = 4.76,
                                   carbon_system_params                          = CarbonSystemParameters(eltype(grid)),
                                   atmospheric_pCO₂                              = 280e-6,
				   )
    FT = eltype(grid)
    
    if incident_PAR isa Number
        surface_PAR = incident_PAR            
        incident_PAR = Field{Center,Center,Nothing}(grid)           
        set!(incident_PAR, surface_PAR)            
        fill_halo_regions!(incident_PAR)
    end
    S = typeof(incident_PAR) # surface type
    if atmospheric_pCO₂ isa Number
        surface_pCO₂ = atmospheric_pCO₂            
        atmospheric_pCO₂ = Field{Center,Center,Nothing}(grid)           
        set!(atmospheric_pCO₂, surface_pCO₂)            
        fill_halo_regions!(atmospheric_pCO₂)
    end

    no_penetration = ImpenetrableBoundaryCondition()
    sinking_bcs = FieldBoundaryConditions(grid, (Center(), Center(), Face()); top=no_penetration, bottom=no_penetration)

    if particulate_organic_phosphorus_sinking_velocity  isa Number
        w₀ᵖᵒᶜ = particulate_organic_phosphorus_sinking_velocity 
        particulate_organic_phosphorus_sinking_velocity  = ZFaceField(grid, boundary_conditions = sinking_bcs)
        set!(particulate_organic_phosphorus_sinking_velocity , w₀ᵖᵒᶜ)
        fill_halo_regions!(particulate_organic_phosphorus_sinking_velocity )
    end

    if particulate_inorganic_carbon_sinking_velocity isa Number
        w₀ᵖⁱᶜ = particulate_inorganic_carbon_sinking_velocity
        particulate_inorganic_carbon_sinking_velocity = ZFaceField(grid, boundary_conditions = sinking_bcs)
        set!(particulate_inorganic_carbon_sinking_velocity, w₀ᵖⁱᶜ)
        fill_halo_regions!(particulate_inorganic_carbon_sinking_velocity)
    end
    W = typeof(particulate_organic_phosphorus_sinking_velocity) # vertical velocity type

    # Set initial first guess of ocean_pH
    ocean_pH = Field{Center,Center,Center}(grid)
    set!(ocean_pH, one(FT)*8)

    # initial guess of calcite saturation state
    calcite_saturation_state = Field{Center,Center,Center}(grid)
    set!(calcite_saturation_state, one(FT))

    # Set initial first guess of ocean_pCO₂
    ocean_pCO₂                 = Field{Center,Center,Center}(grid)
    set!(ocean_pCO₂, atmospheric_pCO₂)
        
    atmospheric_CO₂_solubility = Field{Center,Center,Center}(grid)
    oceanic_CO₂_solubility     = Field{Center,Center,Center}(grid)
    schmidt_number_dic         = Field{Center,Center,Nothing}(grid)
    piston_velocity            = Field{Center,Center,Nothing}(grid)
    surface_wind_speed         = Field{Center,Center,Nothing}(grid)
    CO₂_flux                   = Field{Center,Center,Nothing}(grid)

    PAR                        = Field{Center,Center,Center}(grid)
    NCP                        = Field{Center,Center,Center}(grid)

    # Initialize auxiliary fields to prevent NaN propagation
    set!(atmospheric_CO₂_solubility, one(FT))  # typical seawater CO₂ solubility
    set!(oceanic_CO₂_solubility,     one(FT))
    set!(CO₂_flux,                   zero(FT))
    set!(PAR,                        zero(FT))
    set!(NCP,                        zero(FT))

    fill_halo_regions!(ocean_pH)
    fill_halo_regions!(ocean_pCO₂)
    fill_halo_regions!(calcite_saturation_state)
    fill_halo_regions!(atmospheric_CO₂_solubility)
    fill_halo_regions!(oceanic_CO₂_solubility)
    fill_halo_regions!(schmidt_number_dic)
    fill_halo_regions!(piston_velocity)
    fill_halo_regions!(surface_wind_speed)
    fill_halo_regions!(CO₂_flux)
    fill_halo_regions!(PAR)
    fill_halo_regions!(NCP)

    D = typeof(PAR) # deep ocean/full depth type

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
                                     particulate_inorganic_carbon_sinking_velocity,
                                     convert(FT,calcite_dissolution_rate_saturation_threshold),
                                     convert(FT,calcite_dissolution_rate_above_saturation),    
                                     convert(FT,calcite_dissolution_rate_below_saturation),    
                                     convert(FT,calcite_dissolution_exponent_above_saturation),
                                     convert(FT,calcite_dissolution_exponent_below_saturation),
                                     calcite_saturation_state  ,
                                     PAR                       ,
                                     NCP                       ,
                                     ocean_pH                  ,
                                     ocean_pCO₂                ,
                                     atmospheric_CO₂_solubility,
                                     oceanic_CO₂_solubility    ,
                                     carbon_system_params,
                                     schmidt_number_dic,
                                     piston_velocity,
                                     surface_wind_speed,
                                     atmospheric_pCO₂,
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
                            adapt(to, bgc.particulate_inorganic_carbon_sinking_velocity),
                            adapt(to, bgc.calcite_dissolution_rate_saturation_threshold),
                            adapt(to, bgc.calcite_dissolution_rate_above_saturation),    
                            adapt(to, bgc.calcite_dissolution_rate_below_saturation),    
                            adapt(to, bgc.calcite_dissolution_exponent_above_saturation),
                            adapt(to, bgc.calcite_dissolution_exponent_below_saturation),
                            adapt(to, bgc.calcite_saturation_state  ),
                            adapt(to, bgc.PAR                       ),
                            adapt(to, bgc.NCP                       ),
                            adapt(to, bgc.ocean_pH                  ),
                            adapt(to, bgc.ocean_pCO₂                ),
                            adapt(to, bgc.atmospheric_CO₂_solubility),
                            adapt(to, bgc.oceanic_CO₂_solubility    ),
                            adapt(to, bgc.carbon_system_params),
                            adapt(to, bgc.schmidt_number_dic),
                            adapt(to, bgc.piston_velocity),
                            adapt(to, bgc.surface_wind_speed),
                            adapt(to, bgc.atmospheric_pCO₂),
                            adapt(to, bgc.CO₂_flux),
)

# Define function methods required by Oceananigans BGC interface
const CAN = CarbonAlkalinityNutrients

"""
Required biogeochemical tracers for the CarbonAlkalinityNutrients model
    DIC -> Dissolved Inorganic Carbon
    ALK -> Alkalinity
    PO₄ -> Phosphate (macronutrient)
    NO₃ -> Nitrate (macronutrient)
    DOP -> Dissolved Organic Phosphorus
    POP -> Particulate Organic Phosphorus
    PIC -> Particulate Inorganic Carbon
    Fe  -> Dissolved Iron (micronutrient)
"""
@inline required_biogeochemical_tracers(::CAN) = (:DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :PIC, :Fe)

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
Add a vertical sinking "drift velocity" for the particulate inorganic carbon (PIC) tracer.
"""
@inline function biogeochemical_drift_velocity(bgc::CAN, ::Val{:PIC})
    u = ZeroField()
    v = ZeroField()
    w = bgc.particulate_inorganic_carbon_sinking_velocity
    return (; u, v, w)
end
#@inline biogeochemical_drift_velocity(bgc::CAN, ::Val{:POP})   = bgc.particulate_organic_phosphorus_sinking_velocity
##@inline biogeochemical_advection_scheme(bgc::CAN, ::Val{:POP}) = bgc.advection_scheme

"""
Required biogeochemical auxiliary tracers for the CarbonAlkalinityNutrients model
    PAR: Photosynthetically Active Radiation (W m⁻²)
    ocean_pH: ocean ocean_pH
    ocean_pCO₂: ocean partial pressure of CO₂ (atm)
    atmospheric_CO₂_solubility: solubility of CO₂ in seawater (mol m⁻³ atm⁻¹)
    oceanic_CO₂_solubility: solubility of CO₂ in seawater (mol m⁻³ atm⁻¹)
    CO₂_flux: air-sea flux of CO₂ (mol m⁻² s⁻¹)
"""
@inline required_biogeochemical_auxiliary_fields(::CAN) = tuple(
    :PAR, 
    :NCP,
    :ocean_pH, 
    :calcite_saturation_state,
    :ocean_pCO₂, 
    :atmospheric_CO₂_solubility,
    :oceanic_CO₂_solubility,
    :CO₂_flux,
    :schmidt_number_dic,
    :piston_velocity, 
    :surface_wind_speed,
)
"""
Set the biogeochemical auxiliary tracers for the CarbonAlkalinityNutrients model
"""
@inline biogeochemical_auxiliary_fields(bgc::CAN) = (;
    PAR                        = bgc.PAR,
    NCP                        = bgc.NCP,
    ocean_pH                   = bgc.ocean_pH,
    calcite_saturation_state   = bgc.calcite_saturation_state,
    ocean_pCO₂                 = bgc.ocean_pCO₂, 
    atmospheric_CO₂_solubility = bgc.atmospheric_CO₂_solubility,
    atmospheric_pCO₂           = bgc.atmospheric_pCO₂,
    oceanic_CO₂_solubility     = bgc.oceanic_CO₂_solubility,
    CO₂_flux                   = bgc.CO₂_flux,
    schmidt_number_dic         = bgc.schmidt_number_dic,
    piston_velocity            = bgc.piston_velocity,
    surface_wind_speed         = bgc.surface_wind_speed,
)

## Functions that fill the biogeochemical auxiliary tracers for the CarbonAlkalinityNutrients model
"""
Set the PAR to 40% of the downwelling shortwave radiation
    at the ocean surface.

Use salt forcing to inform DIC and ALK forcing using surface average values
    to convert salinity flux to carbon and alkalinity 'virtual' fluxes.
"""
@kernel function couple_physical_and_bgc_forcing!(
    grid, 
    # Light fields
    sea_ice_concentration, par_fraction, incident_PAR, incoming_shortwave_radiation, 
    # Freshwater flux fields
    surface_DIC, surface_ALK, surface_salinity, salt_forcing, DIC_forcing, ALK_forcing,
)
    FT = eltype(grid)
    i, j = @index(Global, NTuple)
    ks = size(grid, 3)

    # Copy shortwave radiation to incident PAR field at the surface, and scale by par_fraction
    @inbounds if inactive_cell(i, j, ks, grid)
        incident_PAR[i, j, 1] = zero(FT)
    else
        incident_PAR[i, j, 1] = (one(FT) - sea_ice_concentration[i, j, 1]) * par_fraction * incoming_shortwave_radiation[i, j, 1]
    end
    
    # Copy surface salinity flux to DIC and ALK fluxes at the surface, scaled by the surface salinity
    # Does the surface salinity flux already include the effect of sea ice?.
    @inbounds if inactive_cell(i, j, ks, grid) | (surface_salinity[i, j, ks] <= zero(FT))
        DIC_forcing[i, j, 1] = zero(FT)
        ALK_forcing[i, j, 1] = zero(FT)
    else
        Ssafe = 1/max(surface_salinity[i, j, ks], FT(1e-12))
        DIC_forcing[i, j, 1] = surface_DIC[i, j, ks] * Ssafe * salt_forcing[i, j, 1]
        ALK_forcing[i, j, 1] = surface_ALK[i, j, ks] * Ssafe * salt_forcing[i, j, 1]
    end
end

"""
Update Photosynthetically Active Ratiation (PAR) field (i.e. light) for phytoplankton growth
    below the surface based on incident surface PAR and attenuation with depth.
"""
@kernel function update_photosynthetic_active_radiation!(bgc, PAR⁰, PAR, grid) 
    FT = eltype(grid)
    i, j, k = @index(Global, NTuple)

    λ  = bgc.PAR_attenuation_scale
    zᶜ = znode(i, j, k, grid, Center(), Center(), Center())
    inactive = inactive_cell(i, j, k, grid)
    @inbounds PAR[i, j, k] =  ifelse(
		    inactive,
		    zero(FT),
    	    PAR⁰[i, j, 1] * exp(zᶜ / λ),
		    )
end 

"""
    Transfer the surface shortwave radiation and surface salinity forcing, to the incident PAR field
        (par is ~40% of Qs) and surface DIC and ALK forcing fields respectively.
"""
@inline function update_bgc_forcing_fields!(simulation::Simulation)
    ocean = simulation.model.ocean.model
    bgc   = ocean.biogeochemistry
    grid  = ocean.grid

    launch!(
        architecture(grid), grid, :xy, couple_physical_and_bgc_forcing!,
        grid,
        # Light fields
        simulation.model.interfaces.exchanger.sea_ice.state.ℵ,
        bgc.PAR_fraction_of_incoming_solar_radiation,
        bgc.incident_PAR, 
        simulation.model.interfaces.atmosphere_ocean_interface.fluxes.downwelling_shortwave,
        # Freshwater flux fields
        ocean.tracers.DIC,
        ocean.tracers.ALK,
        ocean.tracers.S,
        simulation.model.interfaces.net_fluxes.ocean.S,
        ocean.tracers.DIC.boundary_conditions.top.condition,
        ocean.tracers.ALK.boundary_conditions.top.condition,
    )

    return nothing
end

"""
    solve_ocean_carbon_system!(
            grid,
            solver_params,
            reference_density,
            temperature, 
            salinity, 
            DIC, 
            ALK, 
            PO4, 
            ocean_pH, 
            atmospheric_pCO₂,
            ocean_pCO₂,
            atmospheric_CO₂_solubility, 
            oceanic_CO₂_solubility,
            calcite_saturation_state,
    )

Compute the oceanic pCO₂ using the UniversalRobustCarbonSystem solver.

Arguments:
- `grid::RegularCartesianGrid`: The model grid.
- `solver_params::NamedTuple`: The parameters for the UniversalRobustCarbonSystem solver.
- `reference_density::Float64`: The reference density of seawater.
- `Θᶜ::Field{Center, Center, Nothing}`: Temperature in degrees Celsius.
- `Sᴬ::Field{Center, Center, Nothing}`: Salinity in PSU.
- `Cᵀ::Field{Center, Center, Nothing}`: Total dissolved inorganic carbon (DIC) in mol kg⁻¹.
- `Aᵀ::Field{Center, Center, Nothing}`: Total alkalinity (ALK) in mol kg⁻¹.
- `Pᵀ::Field{Center, Center, Nothing}`: Phosphate concentration in mol kg⁻¹.
- `ocean_pH::Field{Center, Center, Center}`: The computed ocean_pH.
- `atmospheric_CO₂::Field{Center, Center, Nothing}`: The partial pressure of CO₂ in the atmosphere.
- `ocean_pCO₂::Field{Center, Center, Nothing}`: The computed oceanic pCO₂.
- `atmospheric_CO₂_solubility::Field{Center, Center, Nothing}`: The solubility of CO₂ in the atmosphere.
- `oceanic_CO₂_solubility::Field{Center, Center, Nothing}`: The solubility of CO₂ in the ocean.
- `calcite_saturation_state::Field{Center, Center, Center}`: The computed calcite saturation state.
"""
@kernel function solve_ocean_carbon_system!(
    grid,
    solver_params,
    reference_density,
    temperature, 
    salinity, 
    DIC, 
    ALK, 
    PO4, 
    ocean_pH, 
    atmospheric_pCO₂,
    ocean_pCO₂,
    atmospheric_CO₂_solubility, 
    oceanic_CO₂_solubility,
    calcite_saturation_state,
    )

    FT = eltype(grid)
    i, j, k = @index(Global, NTuple)

    if !inactive_cell(i, j, k, grid)
        @inbounds begin
            z      = znode(i, j, k, grid, Center(), Center(), Center())
            Δp_bar = max(zero(FT), -z / FT(10) - FT(1)) # Applied pressure in bar, with 1 bar at the surface
            
	    ## compute oceanic pCO₂ using the UniversalRobustCarbonSystem solver
            CarbonSolved = UniversalRobustCarbonSystem(;
                pH      = ocean_pH[i, j, k],
                pCO₂ᵃᵗᵐ = atmospheric_pCO₂[i, j, 1],
                Θᶜ      = temperature[i, j, k],
                Sᴬ      = max(FT(1e-12), salinity[i, j, k]),  # guard against S ≤ 0 (ice melt)
                Δpᵦₐᵣ   = Δp_bar,
                Cᵀ      = DIC[i, j, k]/reference_density,    # mol/m3 to mol/kg
                Aᵀ      = ALK[i, j, k]/reference_density,    # mol/m3 to mol/kg
                Pᵀ      = PO4[i, j, k]/reference_density,    # mol/m3 to mol/kg
                Siᵀ     = PO4[i, j, k]*15/reference_density, # mol/m3 to mol/kg
                params  = solver_params,
                )

            ocean_pCO₂[i, j, k]                 = ifelse(
                isnan(CarbonSolved.pCO₂ᵒᶜᵉ), 
                ocean_pCO₂[i, j, k],                 
                CarbonSolved.pCO₂ᵒᶜᵉ,
            )
            atmospheric_CO₂_solubility[i, j, k] = ifelse(
                isnan(CarbonSolved.Pᵈⁱᶜₖₛₒₗₐ), 
                atmospheric_CO₂_solubility[i, j, k],
                CarbonSolved.Pᵈⁱᶜₖₛₒₗₐ,
            )
            oceanic_CO₂_solubility[i, j, k]     = ifelse(
                isnan(CarbonSolved.Pᵈⁱᶜₖ₀),
                oceanic_CO₂_solubility[i, j, k],
                CarbonSolved.Pᵈⁱᶜₖ₀,
            )
            ocean_pH[i, j, k]                   = ifelse(
                isnan(CarbonSolved.pH),
                ocean_pH[i, j, k],
                CarbonSolved.pH,
            )
            calcite_saturation_state[i, j, k] = ifelse(
                isnan(CarbonSolved.ΩCᶜᵃˡ),
                calcite_saturation_state[i, j, k],
                CarbonSolved.ΩCᶜᵃˡ,
            )
        end
    else
        @inbounds begin
            ocean_pCO₂[i, j, k]                 = zero(FT)
            atmospheric_CO₂_solubility[i, j, k] = zero(FT)
            oceanic_CO₂_solubility[i, j, k]     = zero(FT)
            ocean_pH[i, j, k]                   = zero(FT)
            calcite_saturation_state[i, j, k]   = zero(FT)
        end
    end
end

"""
    Monod function for nutrient limitation, used in the calculation of net community production (NCP) based on light and nutrient availability.
"""
@inline monod(conc, half_saturationk) = (conc⁺ = max(conc, zero(conc)); conc⁺ / (conc⁺ + half_saturationk))

@kernel function update_net_community_production!(bgc, NCP, grid, PO₄, NO₃, Fe)
    FT = eltype(grid)
    i, j, k = @index(Global, NTuple)

    μᵖ = bgc.maximum_net_community_production_rate
    kᴵ = bgc.PAR_half_saturation
    kᴾ = bgc.phosphate_half_saturation
    kᴺ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation

    @inbounds begin
        inactive = inactive_cell(i, j, k, grid)
        Lₗᵢₘ = monod(bgc.PAR[i, j, k], kᴵ)
        Pₗᵢₘ = monod(PO₄[i, j, k], kᴾ)
        Nₗᵢₘ = monod(NO₃[i, j, k], kᴺ)
        Fₗᵢₘ = monod(Fe[i, j, k], kᶠ)

        NCP[i, j, k] = ifelse(
            inactive,
            zero(FT),
            max(zero(FT), μᵖ * Lₗᵢₘ * min(Pₗᵢₘ, Nₗᵢₘ, Fₗᵢₘ)),
            )
    end
end

@kernel function update_par_and_ncp!(bgc, incident_PAR, PAR, NCP, grid, Δt, PO₄, NO₃, Fe)
    FT = eltype(grid)
    i, j, k = @index(Global, NTuple)

    μᵖ = bgc.maximum_net_community_production_rate
    kᴵ = bgc.PAR_half_saturation
    kᴾ = bgc.phosphate_half_saturation
    kᴺ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    λ  = bgc.PAR_attenuation_scale
    Rᴺ = bgc.stoichoimetric_ratio_nitrate_to_phosphate
    Rᶠ = bgc.stoichoimetric_ratio_iron_to_phosphate

    @inbounds begin
        zᶜ = znode(i, j, k, grid, Center(), Center(), Center())
        inactive = inactive_cell(i, j, k, grid)

        PAR[i, j, k] = ifelse(
            inactive,
            zero(FT),
            incident_PAR[i, j, 1] * exp(zᶜ / λ),
        )

        Lₗᵢₘ = monod(bgc.PAR[i, j, k], kᴵ)
        Pₗᵢₘ = monod(PO₄[i, j, k], kᴾ)
        Nₗᵢₘ = monod(NO₃[i, j, k], kᴺ)
        Fₗᵢₘ = monod(Fe[i, j, k], kᶠ)

        #NCP[i, j, k] = ifelse(
        #    inactive,
        #    zero(FT),
        #    max(zero(FT), μᵖ * Lₗᵢₘ * min(Pₗᵢₘ, Nₗᵢₘ, Fₗᵢₘ)),
        #    )
        # no element may be drawn below zero in one step (γ ≈ 0.5 for AB2 headroom)
        # This is the discrete analogue of the modified-Patankar idea from Burchard, Deleersnijder & Meister (2003), 
        #  specialised to a Liebig-limited uptake = unconditional positivity and exact conservation without an implicit solve.        
        #γ = convert(FT, 0.5)
        γ = one(FT)
        
        ncp_max = γ / Δt * min(max(PO₄[i, j, k], zero(FT)),
                               max(NO₃[i, j, k], zero(FT)) / Rᴺ,
                               max(Fe[ i, j, k], zero(FT)) / Rᶠ)
    
        ncp_local = max(zero(FT), μᵖ * Lₗᵢₘ * min(Pₗᵢₘ, Nₗᵢₘ, Fₗᵢₘ))
        
        NCP[i,j,k] = ifelse(
            inactive,
            zero(FT),
            min(ncp_local, ncp_max),
        )
    end
end

"""
Update BGC auxiliary fields
"""
@inline function update_biogeochemical_state!(bgc::CAN, model)
    arch = architecture(model.grid)
    grid = model.grid
    FT = eltype(grid)

    # Multi-stage steppers call update_state! multiple times per timestep.
    # Avoid redundant carbon-system solves by running it once per full step.
    stage_count = if hasproperty(model.timestepper, :β)
        length(getproperty(model.timestepper, :β))
    elseif hasproperty(model.timestepper, :γ³)
        3
    else
        0
    end
    if (model.clock.iteration == 0) || (stage_count == 0) || (model.clock.stage == stage_count)
        # Fused PAR+NCP update to reduce memory traffic and kernel launches.
        kernel_args = (
            bgc,
            bgc.incident_PAR,
            bgc.PAR,
            bgc.NCP,
            grid,
            model.clock.last_Δt,
            model.tracers.PO₄, 
            model.tracers.NO₃, 
            model.tracers.Fe,
        )

        launch!(arch, 
                grid, 
                :xyz, 
                update_par_and_ncp!,
                kernel_args...,
        )

        # Solve the ocean carbon system and calculate ocean_pH, pCO₂, etc.
        # This is intentionally run once per full timestep for multi-stage schemes.
        kernel_args = (
            grid,
            bgc.carbon_system_params,
            bgc.reference_density,
            model.tracers.T,
            model.tracers.S,
            model.tracers.DIC,
            model.tracers.ALK,
            model.tracers.PO₄,
            bgc.ocean_pH,
            bgc.atmospheric_pCO₂,
            bgc.ocean_pCO₂,
            bgc.atmospheric_CO₂_solubility, 
            bgc.oceanic_CO₂_solubility,
            bgc.calcite_saturation_state,
        )

        launch!(
            arch,
            grid,
            :xyz,
            solve_ocean_carbon_system!,
            kernel_args...,
        )
    end
    return nothing
end

"""
Update BGC tendencies (after everything is calculated) for the CarbonAlkalinityNutrients model.
"""
@inline function update_biogeochemical_tendencies!(bgc::CAN, model)
    arch = architecture(model.grid)
    grid = model.grid

    return nothing
end

## Now define functions that relate to the specific BGC model
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
    dissolved_organic_phosphorus_remin(grid, remineralization_rate, 
                                      dissolved_organic_phosphorus_concentration)
Calculate the remineralization of dissolved organic phosphorus.
"""
@inline function dissolved_organic_phosphorus_remin(grid, 
                                                    remineralization_rate, 
                                                    dissolved_organic_phosphorus_concentration) 
        FT = eltype(grid)
    return max(zero(FT), remineralization_rate * dissolved_organic_phosphorus_concentration)
end

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
        FT = eltype(grid)
        Rᵣ = option_of_particulate_remin                                  
        r = particulate_organic_phosphorus_remin_timescale
        b = martin_curve_exponent    
        wₛ = particulate_organic_phosphorus_sinking_velocity
        λ = PAR_attenuation_scale
        z  = depth
        fᵢ= percent_light
        z₀ = log(fᵢ)*λ # The base of the euphotic layer depth (z₀) where PAR is degraded down to 1%     
        POP = particulate_organic_phosphorus_concentration

        γ = ifelse(Rᵣ == one(FT), b * wₛ / (z + z₀), r)
    return γ * POP
#    return ifelse(Rᵣ == one(FT), max(zero(FT), b * wₛ / (z + z₀) * POP), max(zero(FT), r * POP))
end


"""
Calculate remineralization of particulate inorganic carbon using Naviaux et al. 2019,
 Marine Chemistry dissolution rates. The value of 0.8272 differs slightly from the paper (0.8),
 but reduces the large discontinuity in dissolution rates between the two states.

Default dissolution rate  values are slightly different from the paper, but consistent with
 the Darwin model. Values include a ~1e5 cm2/mol conversion factor for specific surface area (SSA) 
 of CaCO3 (e.g. benthic forams or synthetic calcite), for example, Subhas et al. 2018 (Marine Chemistry).
 dissolution_rate = tuple(5.22e-9, 0.11, 1.65e-5, 4.76)
"""
@inline function particulate_inorganic_carbon_remin(
                            grid, 
                            calcite_dissolution_rate_saturation_threshold,
                            calcite_dissolution_rate_above_saturation,    
                            calcite_dissolution_rate_below_saturation,    
                            calcite_dissolution_exponent_above_saturation,
                            calcite_dissolution_exponent_below_saturation,
                            PIC, 
                            ΩC)
    FT = eltype(grid)
    ΩCᵗʰʳᵉˢʰ = calcite_dissolution_rate_saturation_threshold
    ΩC⁺      = calcite_dissolution_rate_above_saturation
    ΩC⁻      = calcite_dissolution_rate_below_saturation
    ΩC⁺ᵉˣᵖ   = calcite_dissolution_exponent_above_saturation
    ΩC⁻ᵉˣᵖ   = calcite_dissolution_exponent_below_saturation

    γ = ifelse(
        ΩC > ΩCᵗʰʳᵉˢʰ,
        max(zero(FT), ΩC⁺ * max(zero(FT), one(FT) - ΩC)^ΩC⁺ᵉˣᵖ),
        max(zero(FT), ΩC⁻ * max(zero(FT), one(FT) - ΩC)^ΩC⁻ᵉˣᵖ),
    )
    return γ * PIC
#    return ifelse(
#        ΩC > ΩCᵗʰʳᵉˢʰ,
#        max(zero(FT), PIC * ΩC⁺ * max(zero(FT), one(FT) - ΩC)^ΩC⁺ᵉˣᵖ),
#        max(zero(FT), PIC * ΩC⁻ * max(zero(FT), one(FT) - ΩC)^ΩC⁻ᵉˣᵖ),
#    )
end

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
    λ     = bgc.PAR_attenuation_scale
    fᵢ    = bgc.PAR_percent
    γ     = bgc.dissolved_organic_phosphorus_remin_timescale
    r     = bgc.particulate_organic_phosphorus_remin_timescale
    α     = bgc.fraction_of_particulate_export
    Rᶜᴾ   = bgc.stoichoimetric_ratio_carbon_to_phosphate
    Rᶜᵃᶜᵒ³ = bgc.rain_ratio_inorganic_to_organic_carbon
    b     = bgc.martin_curve_exponent
    wₛ    = bgc.particulate_organic_phosphorus_sinking_velocity
    Rᵣ    = bgc.option_of_particulate_remin
    ΩCᵗʰʳᵉˢʰ = bgc.calcite_dissolution_rate_saturation_threshold
    ΩC⁺      = bgc.calcite_dissolution_rate_above_saturation
    ΩC⁻      = bgc.calcite_dissolution_rate_below_saturation
    ΩC⁺ᵉˣᵖ   = bgc.calcite_dissolution_exponent_above_saturation
    ΩC⁻ᵉˣᵖ   = bgc.calcite_dissolution_exponent_below_saturation

    @inbounds begin
       z         = znode(i, j, k, grid, Center(), Center(), Center())
       NCP       = bgc.NCP[i, j, k]
       DOP       = fields.DOP[i, j, k]
       POP       = fields.POP[i, j, k]
       PIC       = fields.PIC[i, j, k]
       ΩC        = bgc.calcite_saturation_state[i, j, k]
       kface_top = min(k + 1, size(wₛ, 3))
       Wₛ        = (wₛ[i, j, k] + wₛ[i, j, kface_top]) / 2
    end

    return Rᶜᴾ * (
                    dissolved_organic_phosphorus_remin(grid, γ, DOP) - 
                    NCP -
                    Rᶜᵃᶜᵒ³ * α * NCP +
                    particulate_organic_phosphorus_remin(grid, Rᵣ, r, b, Wₛ, λ, z, fᵢ, POP)) +
              particulate_inorganic_carbon_remin(grid, ΩCᵗʰʳᵉˢʰ, ΩC⁺, ΩC⁻, ΩC⁺ᵉˣᵖ, ΩC⁻ᵉˣᵖ, PIC, ΩC)
end

"""
Tracer sources and sinks for alkalinity (ALK)
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:ALK}, clock, fields)
    λ      = bgc.PAR_attenuation_scale
    fᵢ     = bgc.PAR_percent
    γ      = bgc.dissolved_organic_phosphorus_remin_timescale
    r      = bgc.particulate_organic_phosphorus_remin_timescale
    α      = bgc.fraction_of_particulate_export
    Rᴺᴾ    = bgc.stoichoimetric_ratio_nitrate_to_phosphate
    Rᶜᴾ    = bgc.stoichoimetric_ratio_carbon_to_phosphate
    Rᶜᵃᶜᵒ³ = bgc.rain_ratio_inorganic_to_organic_carbon
    b      = bgc.martin_curve_exponent
    wₛ     = bgc.particulate_organic_phosphorus_sinking_velocity
    Rᵣ     = bgc.option_of_particulate_remin
    ΩCᵗʰʳᵉˢʰ = bgc.calcite_dissolution_rate_saturation_threshold
    ΩC⁺      = bgc.calcite_dissolution_rate_above_saturation
    ΩC⁻      = bgc.calcite_dissolution_rate_below_saturation
    ΩC⁺ᵉˣᵖ   = bgc.calcite_dissolution_exponent_above_saturation
    ΩC⁻ᵉˣᵖ   = bgc.calcite_dissolution_exponent_below_saturation

    @inbounds begin
       z         = znode(i, j, k, grid, Center(), Center(), Center())
       NCP       = bgc.NCP[i, j, k]
       DOP       = fields.DOP[i, j, k]
       POP       = fields.POP[i, j, k]
       PIC       = fields.PIC[i, j, k]
       ΩC        = bgc.calcite_saturation_state[i, j, k]
       kface_top = min(k + 1, size(wₛ, 3))
       Wₛ        = (wₛ[i, j, k] + wₛ[i, j, kface_top]) / 2
    end

    # Nitrate uptake/remineralization releases/consumes alkalinity via Rᴺᴾ,
    # while CaCO₃ formation/dissolution consumes/releases 2 mol alkalinity per
    # mol of inorganic carbon (Rᶜᴾ * Rᶜᵃᶜᵒ³ * α * NCP).
    return 2 * ( 
                    particulate_inorganic_carbon_remin(grid, ΩCᵗʰʳᵉˢʰ, ΩC⁺, ΩC⁻, ΩC⁺ᵉˣᵖ, ΩC⁻ᵉˣᵖ, PIC, ΩC) - 
                    Rᶜᴾ * Rᶜᵃᶜᵒ³ * α * NCP
           ) -Rᴺᴾ * (
                       dissolved_organic_phosphorus_remin(grid, γ, DOP) +
                       particulate_organic_phosphorus_remin(grid, Rᵣ, r, b, Wₛ, λ, z, fᵢ, POP) -
                       NCP
           )
end

"""
Tracer sources and sinks for dissolved inorganic phosphate (PO₄)
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:PO₄}, clock, fields)
    λ  = bgc.PAR_attenuation_scale
    fᵢ = bgc.PAR_percent
    γ  = bgc.dissolved_organic_phosphorus_remin_timescale
    r  = bgc.particulate_organic_phosphorus_remin_timescale
    b  = bgc.martin_curve_exponent
    wₛ = bgc.particulate_organic_phosphorus_sinking_velocity
    Rᵣ = bgc.option_of_particulate_remin

    @inbounds begin
       z         = znode(i, j, k, grid, Center(), Center(), Center())
       NCP       = bgc.NCP[i, j, k]
       DOP       = fields.DOP[i, j, k]
       POP       = fields.POP[i, j, k]
       kface_top = min(k + 1, size(wₛ, 3))
       Wₛ        = (wₛ[i, j, k] + wₛ[i, j, kface_top]) / 2
    end

    return dissolved_organic_phosphorus_remin(grid, γ, DOP) +
           particulate_organic_phosphorus_remin(grid, Rᵣ, r, b, Wₛ, λ, z, fᵢ, POP) -
           NCP
end

"""
Tracer sources and sinks for dissolved inorganic nitrate (NO₃)
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:NO₃}, clock, fields)
    λ   = bgc.PAR_attenuation_scale
    fᵢ  = bgc.PAR_percent
    γ   = bgc.dissolved_organic_phosphorus_remin_timescale
    r   = bgc.particulate_organic_phosphorus_remin_timescale
    Rᴺᴾ = bgc.stoichoimetric_ratio_nitrate_to_phosphate
    b   = bgc.martin_curve_exponent
    wₛ  = bgc.particulate_organic_phosphorus_sinking_velocity
    Rᵣ  = bgc.option_of_particulate_remin

    @inbounds begin
       z         = znode(i, j, k, grid, Center(), Center(), Center())
       NCP       = bgc.NCP[i, j, k]
       DOP       = fields.DOP[i, j, k]
       POP       = fields.POP[i, j, k]
       kface_top = min(k + 1, size(wₛ, 3))
       Wₛ        = (wₛ[i, j, k] + wₛ[i, j, kface_top]) / 2
    end

    return Rᴺᴾ * (
           dissolved_organic_phosphorus_remin(grid, γ, DOP) +
           particulate_organic_phosphorus_remin(grid, Rᵣ, r, b, Wₛ, λ, z, fᵢ, POP) -
           NCP
           )
end

"""
Tracer sources and sinks for dissolved iron (FeT)
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:Fe}, clock, fields)
    λ     = bgc.PAR_attenuation_scale
    fᵢ    = bgc.PAR_percent
    γ     = bgc.dissolved_organic_phosphorus_remin_timescale
    r     = bgc.particulate_organic_phosphorus_remin_timescale
    Rᶠᴾ   = bgc.stoichoimetric_ratio_iron_to_phosphate
    Lᶠᵉ   = bgc.ligand_concentration
    β     = bgc.ligand_stability_coefficient
    kˢᶜᵃᵛ = bgc.iron_scavenging_rate
    b     = bgc.martin_curve_exponent
    wₛ    = bgc.particulate_organic_phosphorus_sinking_velocity
    Rᵣ    = bgc.option_of_particulate_remin

    @inbounds begin
       z         = znode(i, j, k, grid, Center(), Center(), Center())
       NCP       = bgc.NCP[i, j, k]
       Feₜ       = fields.Fe[i, j, k]
       DOP       = fields.DOP[i, j, k]
       POP       = fields.POP[i, j, k]
       kface_top = min(k + 1, size(wₛ, 3))
       Wₛ        = (wₛ[i, j, k] + wₛ[i, j, kface_top]) / 2
    end

    return Rᶠᴾ * (
                  dissolved_organic_phosphorus_remin(grid, γ, DOP) +
                  particulate_organic_phosphorus_remin(grid, Rᵣ, r, b, Wₛ, λ, z, fᵢ, POP) -
                  NCP
                 ) +
                 iron_sources() -
                 iron_scavenging(kˢᶜᵃᵛ, Feₜ, Lᶠᵉ, β)
end


"""
Tracer sources and sinks for dissolved organic phosphorus (DOP)
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:DOP}, clock, fields)
    γ = bgc.dissolved_organic_phosphorus_remin_timescale
    α = bgc.fraction_of_particulate_export

    @inbounds begin
       NCP = bgc.NCP[i, j, k]
       DOP = fields.DOP[i, j, k]
    end

    return (1 - α) * NCP - dissolved_organic_phosphorus_remin(grid, γ, DOP)
end

"""
Tracer sources and sinks for Particulate Organic Phosphorus (POP).
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:POP}, clock, fields)
    λ  = bgc.PAR_attenuation_scale
    fᵢ = bgc.PAR_percent
    r  = bgc.particulate_organic_phosphorus_remin_timescale
    α  = bgc.fraction_of_particulate_export
    b  = bgc.martin_curve_exponent
    wₛ = bgc.particulate_organic_phosphorus_sinking_velocity
    Rᵣ = bgc.option_of_particulate_remin

    @inbounds begin
       z         = znode(i, j, k, grid, Center(), Center(), Center())
       NCP       = bgc.NCP[i, j, k]
       POP       = fields.POP[i, j, k]
       kface_top = min(k + 1, size(wₛ, 3))
       Wₛ        = (wₛ[i, j, k] + wₛ[i, j, kface_top]) / 2
    end

    return α * NCP - particulate_organic_phosphorus_remin(grid, Rᵣ, r, b, Wₛ, λ, z, fᵢ, POP)
end

"""
Tracer sources and sinks for Particulate Inorganic Carbon (PIC).
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:PIC}, clock, fields)
    α      = bgc.fraction_of_particulate_export
    Rᶜᴾ    = bgc.stoichoimetric_ratio_carbon_to_phosphate
    Rᶜᵃᶜᵒ³ = bgc.rain_ratio_inorganic_to_organic_carbon
    ΩCᵗʰʳᵉˢʰ = bgc.calcite_dissolution_rate_saturation_threshold
    ΩC⁺      = bgc.calcite_dissolution_rate_above_saturation
    ΩC⁻      = bgc.calcite_dissolution_rate_below_saturation
    ΩC⁺ᵉˣᵖ   = bgc.calcite_dissolution_exponent_above_saturation
    ΩC⁻ᵉˣᵖ   = bgc.calcite_dissolution_exponent_below_saturation

    @inbounds begin
        NCP    = bgc.NCP[i, j, k]
        PIC    = fields.PIC[i, j, k]
        ΩC     = bgc.calcite_saturation_state[i, j, k]
    end

    return Rᶜᴾ * Rᶜᵃᶜᵒ³ * α * NCP - particulate_inorganic_carbon_remin(grid, ΩCᵗʰʳᵉˢʰ, ΩC⁺, ΩC⁻, ΩC⁺ᵉˣᵖ, ΩC⁻ᵉˣᵖ, PIC, ΩC)
end
