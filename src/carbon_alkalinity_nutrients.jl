using Oceananigans.Units: days
using Oceananigans.Grids: znode, Center
using Oceananigans.Fields: ZeroField, ConstantField
using Oceananigans.Biogeochemistry: AbstractBiogeochemistry

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers, biogeochemical_drift_velocity

const c = Center()

struct CarbonAlkalinityNutrients{FT, W} <: AbstractBiogeochemistry
    reference_density                             :: FT
    maximum_net_community_production_rate         :: FT # mol PO₄ m⁻³ s⁻¹
    phosphate_half_saturation                     :: FT # mol PO₄ m⁻³
    nitrate_half_saturation                       :: FT # mol NO₃ m⁻³
    iron_half_saturation                          :: FT # mol Fe m⁻³
    PAR_half_saturation                           :: FT  # W m⁻²
    PAR_attenuation_scale                         :: FT  # m
    fraction_of_particulate_export                :: FT
    dissolved_organic_phosphate_remin_timescale   :: FT # s⁻¹
    stoichoimetric_ratio_carbon_to_phosphate      :: FT 
    stoichoimetric_ratio_nitrate_to_phosphate     :: FT 
    stoichoimetric_ratio_phosphate_to_oxygen      :: FT 
    stoichoimetric_ratio_phosphate_to_iron        :: FT 
    stoichoimetric_ratio_carbon_to_nitrate        :: FT 
    stoichoimetric_ratio_carbon_to_oxygen         :: FT 
    stoichoimetric_ratio_carbon_to_iron           :: FT 
    stoichoimetric_ratio_silicate_to_phosphate    :: FT
    particulate_organic_phosphate_remin_timescale :: FT 
    rain_ratio_inorganic_to_organic_carbon        :: FT  
    iron_scavenging_rate                          :: FT # s⁻¹
    ligand_concentration                          :: FT # mol L m⁻³
    ligand_stability_coefficient                  :: FT
    particulate_organic_phosphate_sinking_speed   :: W  # m s⁻¹       
end

"""
    CarbonAlkalinityNutrients(  reference_density                             = 1024.5,
                                maximum_net_community_production_rate         = 1 / day,
                                phosphate_half_saturation                     = 1e-7 * reference_density,
                                nitrate_half_saturation                       = 1.6e-6 * reference_density,
                                iron_half_saturation                          = 1e-10 * reference_density,
                                PAR_half_saturation                           = 10.0,
                                PAR_attenuation_scale                         = 25.0,
                                fraction_of_particulate_export                = 0.33
                                dissolved_organic_phosphate_remin_timescale   = 1 / 30day,
                                stoichoimetric_ratio_carbon_to_phosphate      = 106.0
                                stoichoimetric_ratio_nitrate_to_phosphate     = 16.0
                                stoichoimetric_ratio_phosphate_to_oxygen      = 170.0,
                                stoichoimetric_ratio_phosphate_to_iron        = 4.68e-4
                                stoichoimetric_ratio_carbon_to_nitrate        = 106 / 16
                                stoichoimetric_ratio_carbon_to_oxygen         = 106 / 170,
                                stoichoimetric_ratio_carbon_to_iron           = 106 / 1.e-3,
                                stoichoimetric_ratio_silicate_to_phosphate    = 15.0,
                                rain_ratio_inorganic_to_organic_carbon        = 1e-1,
                                particulate_organic_phosphate_remin_timescale = 0.03 / day,
                                iron_scavenging_rate                          = 5e-4 / day,
                                ligand_concentration                          = 1e-9 * reference_density,
                                ligand_stability_coefficient                  = 1e8
                                particulate_organic_phosphate_sinking_speed   = -10.0 / day)

Return a seven-tracer biogeochemistry model for the interaction of carbon, alkalinity, and nutrients.

Keyword Arguments
=================

Tracer names
============
* `DIC`: Dissolved Inorganic Carbon

* `ALK`: Alkalinity

* `PO₄`: Phosphate (macronutrient)

* `NO₃`: Nitrate (macronutrient)

* `DOP`: Dissolved Organic Phosphate (macronutrient)

* `Fe`: Dissolved Iron (micronutrient)

Biogeochemical functions
========================
* transitions for `DIC`, `ALK`, `PO₄`, `NO₃`, `POP`, `DOP`, `POP` and `Fe`

* `biogeochemical_drift_velocity` for `DOP`, modeling the sinking of detritus at
  a constant `detritus_sinking_speed`.
"""
function CarbonAlkalinityNutrients(; grid,
                                   reference_density                            = 1024.5,
                                   maximum_net_community_production_rate        = 2.e-3 / 365.25days, # mol PO₄ m⁻³ s⁻¹
                                   phosphate_half_saturation                    = 5.e-7 * reference_density, # mol PO₄ m⁻³
                                   nitrate_half_saturation                      = 7.e-6 * reference_density, # mol NO₃ m⁻³
                                   iron_half_saturation                         = 1.e-10 * reference_density, # mol Fe m⁻³
                                   PAR_half_saturation                          = 30.0,  # W m⁻²
                                   PAR_attenuation_scale                        = 25.0,  # m
                                   fraction_of_particulate_export               = 0.33,
                                   dissolved_organic_phosphate_remin_timescale  = 2. / 365.25days, # s⁻¹
                                   stoichoimetric_ratio_carbon_to_phosphate     = 117.0,
                                   stoichoimetric_ratio_nitrate_to_phosphate    = 16.0,
                                   stoichoimetric_ratio_phosphate_to_oxygen     = 170.0, 
                                   stoichoimetric_ratio_phosphate_to_iron       = 4.68e-4,
                                   stoichoimetric_ratio_carbon_to_nitrate       = 117. / 16.,
                                   stoichoimetric_ratio_carbon_to_oxygen        = 117. / 170., 
                                   stoichoimetric_ratio_carbon_to_iron          = 117. / 4.68e-4,
                                   stoichoimetric_ratio_silicate_to_phosphate   = 15.0,
                                   rain_ratio_inorganic_to_organic_carbon       = 1.e-2,
                                   particulate_organic_phosphate_remin_timescale= 0.03 / day, 
                                   iron_scavenging_rate                         = 0.2 / 365.25days, # s⁻¹
                                   ligand_concentration                         = 1e-9 * reference_density, # mol L m⁻³
                                   ligand_stability_coefficient                 = 1e8,
                                   particulate_organic_phosphate_sinking_speed  = -10.0 / day)

    if particulate_organic_phosphate_sinking_speed isa Number
        w₀ = particulate_organic_phosphate_sinking_speed
        no_penetration = ImpenetrableBoundaryCondition()

        bcs = FieldBoundaryConditions(grid, (Center, Center, Face),
                                      top=no_penetration, bottom=no_penetration)

        particulate_organic_phosphate_sinking_speed = ZFaceField(grid, boundary_conditions = bcs)

        set!(particulate_organic_phosphate_sinking_speed, w₀)

        fill_halo_regions!(particulate_organic_phosphate_sinking_speed)
    end
                            
    FT = eltype(grid)

    return CarbonAlkalinityNutrients(convert(FT, reference_density),
                                     convert(FT, maximum_net_community_production_rate),
                                     convert(FT, phosphate_half_saturation),
                                     convert(FT, nitrate_half_saturation),
                                     convert(FT, iron_half_saturation),
                                     convert(FT, PAR_half_saturation),
                                     convert(FT, PAR_attenuation_scale),
                                     convert(FT, fraction_of_particulate_export),
                                     convert(FT, dissolved_organic_phosphate_remin_timescale),
                                     convert(FT, stoichoimetric_ratio_carbon_to_phosphate),
                                     convert(FT, stoichoimetric_ratio_nitrate_to_phosphate),
                                     convert(FT, stoichoimetric_ratio_phosphate_to_oxygen),
                                     convert(FT, stoichoimetric_ratio_phosphate_to_iron),
                                     convert(FT, stoichoimetric_ratio_carbon_to_nitrate),
                                     convert(FT, stoichoimetric_ratio_carbon_to_oxygen),
                                     convert(FT, stoichoimetric_ratio_carbon_to_iron),
                                     convert(FT, stoichoimetric_ratio_silicate_to_phosphate),
                                     convert(FT, rain_ratio_inorganic_to_organic_carbon),
                                     convert(FT, particulate_organic_phosphate_remin_timescale),
                                     convert(FT, iron_scavenging_rate),
                                     convert(FT, ligand_concentration),
                                     convert(FT, ligand_stability_coefficient),
                                     particulate_organic_phosphate_sinking_speed)
end

const CAN = CarbonAlkalinityNutrients

@inline required_biogeochemical_tracers(::CAN) = (:DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe)

@inline function biogeochemical_drift_velocity(bgc::CAN, ::Val{:POP})
    u = ZeroField()
    v = ZeroField()
    w = bgc.particulate_organic_phosphate_sinking_speed
    return (; u, v, w)
end

@inline function net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, PO₄, NO₃, Feₜ)
    return μᵖ * I / (I + kᴵ) * min(PO₄ / (PO₄ + kᴾ), NO₃ / (NO₃ + kᴺ), Feₜ / (Feₜ + kᶠ))
end

@inline dissolved_organic_phosphate_remin(γ, DOP) = γ * DOP

@inline particulate_organic_phosphate_remin(r, POP) = r * POP

# exponential remineralization or below the lysocline
@inline particulate_inorganic_carbon_remin() = 0.0

#@inline air_sea_flux_co2() = 0.0

#@inline freshwater_virtual_flux() = 0.0

"""
Iron scavenging should depend on free iron, involves solving a quadratic equation in terms
of ligand concentration and stability coefficient, but this is a simple first order approximation.
"""
@inline function iron_scavenging(kˢᶜᵃᵛ, Fₜ, Lₜ, β)
    # solve for the equilibrium free iron concentration
       # β = FeL / (Feᶠʳᵉᵉ * Lᶠʳᵉᵉ)
       # Lₜ = FeL + Lᶠʳᵉᵉ
       # Fₜ = FeL + Feᶠʳᵉᵉ
       # --> R₁(Feᶠʳᵉᵉ)² + R₂ Feᶠʳᵉᵉ + R₃ = 0
       β⁻¹ = 1/β
       R₁  = 1
       R₂  = (Lₜ + β⁻¹ - Fₜ) 
       R₃  = -(Fₜ * β⁻¹) 

       # simple quadratic solution for roots
       discriminant = ( R₂*R₂ - ( 4*R₁*R₃ ))^(1/2)

       # directly solve for the free iron concentration
       Feᶠʳᵉᵉ = (-R₂ + discriminant) / (2*R₁) 

       # return the linear scavenging rate (net scavenging)
       return (kˢᶜᵃᵛ * Feᶠʳᵉᵉ)
end

@inline iron_sources() = 1e-7

"""
Tracer sources and sinks for DIC
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:DIC}, clock, fields)
    μᵖ = bgc.maximum_net_community_production_rate
    kᴺ = bgc.phosphate_half_saturation
    kᴾ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    kᴵ = bgc.PAR_half_saturation
    λ = bgc.PAR_attenuation_scale
    γ = bgc.dissolved_organic_phosphate_remin_timescale
    r = bgc.particulate_organic_phosphate_remin_timescale
    α = bgc.fraction_of_particulate_export
    Rᶜᴾ = bgc.stoichoimetric_ratio_carbon_to_phosphate   
    Rᴺᴾ = bgc.stoichoimetric_ratio_nitrate_to_phosphate  
    Rᴾᴼ = bgc.stoichoimetric_ratio_phosphate_to_oxygen   
    Rᶜᴺ = bgc.stoichoimetric_ratio_carbon_to_nitrate     
    Rᶜᴼ = bgc.stoichoimetric_ratio_carbon_to_oxygen      
    Rᶜᶠ = bgc.stoichoimetric_ratio_carbon_to_iron        
    Rᶜᵃᶜᵒ³ = bgc.rain_ratio_inorganic_to_organic_carbon     

    # Available photosynthetic radiation
    z = znode(i, j, k, grid, c, c, c)
    # TODO: design a user interface for prescribing incoming shortwave
    I = 700 * exp(z / λ)

    PO₄ = @inbounds fields.PO₄[i, j, k]
    NO₃ = @inbounds fields.NO₃[i, j, k]
    Feₜ = @inbounds fields.Fe[i, j, k]
    DOP = @inbounds fields.DOP[i, j, k]
    POP = @inbounds fields.POP[i, j, k]
    
    return Rᶜᴾ * (
                    dissolved_organic_phosphate_remin(γ, DOP) +
                    particulate_organic_phosphate_remin(r, POP) -
                    (1 + Rᶜᵃᶜᵒ³ * α) * net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, PO₄, NO₃, Feₜ)
                ) +
           particulate_inorganic_carbon_remin()
end

"""
Tracer sources and sinks for ALK
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:ALK}, clock, fields)
    μᵖ = bgc.maximum_net_community_production_rate
    kᴺ = bgc.phosphate_half_saturation
    kᴾ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    kᴵ = bgc.PAR_half_saturation
    λ = bgc.PAR_attenuation_scale
    γ = bgc.dissolved_organic_phosphate_remin_timescale
    r = bgc.dissolved_organic_phosphate_remin_timescale
    α = bgc.fraction_of_particulate_export
    Rᶜᴾ = bgc.stoichoimetric_ratio_carbon_to_phosphate   
    Rᴺᴾ = bgc.stoichoimetric_ratio_nitrate_to_phosphate  
    Rᴾᴼ = bgc.stoichoimetric_ratio_phosphate_to_oxygen   
    Rᶜᴺ = bgc.stoichoimetric_ratio_carbon_to_nitrate     
    Rᶜᴼ = bgc.stoichoimetric_ratio_carbon_to_oxygen      
    Rᶜᶠ = bgc.stoichoimetric_ratio_carbon_to_iron        
    Rᶜᵃᶜᵒ³ = bgc.rain_ratio_inorganic_to_organic_carbon     

    # Available photosynthetic radiation
    z = znode(i, j, k, grid, c, c, c)
    # TODO: design a user interface for prescribing incoming shortwave
    I = 700 * exp(z / λ)

    PO₄ = @inbounds fields.PO₄[i, j, k]
    NO₃ = @inbounds fields.NO₃[i, j, k]
    Feₜ = @inbounds fields.Fe[i, j, k]
    DOP = @inbounds fields.DOP[i, j, k]
    POP = @inbounds fields.POP[i, j, k]
    
    return -Rᴺᴾ * (
        - (1 + Rᶜᵃᶜᵒ³ * α) * net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, PO₄, NO₃, Feₜ) +
        dissolved_organic_phosphate_remin(γ, DOP) +
        particulate_organic_phosphate_remin(r, POP)) +
        2 * particulate_inorganic_carbon_remin()
end

"""
Tracer sources and sinks for PO₄
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:PO₄}, clock, fields)
    μᵖ = bgc.maximum_net_community_production_rate
    kᴺ = bgc.phosphate_half_saturation
    kᴾ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    kᴵ = bgc.PAR_half_saturation
    λ = bgc.PAR_attenuation_scale
    γ = bgc.dissolved_organic_phosphate_remin_timescale
    r = bgc.particulate_organic_phosphate_remin_timescale
    α = bgc.fraction_of_particulate_export
    Rᶜᴾ = bgc.stoichoimetric_ratio_carbon_to_phosphate   
    Rᴺᴾ = bgc.stoichoimetric_ratio_nitrate_to_phosphate  
    Rᴾᴼ = bgc.stoichoimetric_ratio_phosphate_to_oxygen   
    Rᶜᴺ = bgc.stoichoimetric_ratio_carbon_to_nitrate     
    Rᶜᴼ = bgc.stoichoimetric_ratio_carbon_to_oxygen      
    Rᶜᶠ = bgc.stoichoimetric_ratio_carbon_to_iron        
    Rᶜᵃᶜᵒ³ = bgc.rain_ratio_inorganic_to_organic_carbon     

    # Available photosynthetic radiation
    z = znode(i, j, k, grid, c, c, c)
    # TODO: design a user interface for prescribing incoming shortwave
    I = 700 * exp(z / λ)

    PO₄ = @inbounds fields.PO₄[i, j, k]
    NO₃ = @inbounds fields.NO₃[i, j, k]
    Feₜ = @inbounds fields.Fe[i, j, k]
    DOP = @inbounds fields.DOP[i, j, k]
    POP = @inbounds fields.POP[i, j, k]

    return - net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, PO₄, NO₃, Feₜ) +
           dissolved_organic_phosphate_remin(γ, DOP) +
           particulate_organic_phosphate_remin(r, POP)
end

"""
Tracer sources and sinks for NO₃
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:NO₃}, clock, fields)
    μᵖ = bgc.maximum_net_community_production_rate
    kᴺ = bgc.phosphate_half_saturation
    kᴾ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    kᴵ = bgc.PAR_half_saturation
    λ = bgc.PAR_attenuation_scale
    γ = bgc.dissolved_organic_phosphate_remin_timescale
    r = bgc.particulate_organic_phosphate_remin_timescale
    α = bgc.fraction_of_particulate_export
    Rᶜᴾ = bgc.stoichoimetric_ratio_carbon_to_phosphate   
    Rᴺᴾ = bgc.stoichoimetric_ratio_nitrate_to_phosphate  
    Rᴾᴼ = bgc.stoichoimetric_ratio_phosphate_to_oxygen   
    Rᶜᴺ = bgc.stoichoimetric_ratio_carbon_to_nitrate     
    Rᶜᴼ = bgc.stoichoimetric_ratio_carbon_to_oxygen      
    Rᶜᶠ = bgc.stoichoimetric_ratio_carbon_to_iron        
    Rᶜᵃᶜᵒ³ = bgc.rain_ratio_inorganic_to_organic_carbon     

    # Available photosynthetic radiation
    z = znode(i, j, k, grid, c, c, c)
    # TODO: design a user interface for prescribing incoming shortwave
    I = 700 * exp(z / λ)

    PO₄ = @inbounds fields.PO₄[i, j, k]
    NO₃ = @inbounds fields.NO₃[i, j, k]
    Feₜ = @inbounds fields.Fe[i, j, k]
    DOP = @inbounds fields.DOP[i, j, k]
    POP = @inbounds fields.POP[i, j, k]

    return Rᴺᴾ * (
           - net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, PO₄, NO₃, Feₜ) +
           dissolved_organic_phosphate_remin(γ, DOP) +
           particulate_organic_phosphate_remin(r, POP))
end

"""
Tracer sources and sinks for FeT
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:Fe}, clock, fields)
    μᵖ = bgc.maximum_net_community_production_rate
    kᴺ = bgc.phosphate_half_saturation
    kᴾ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    kᴵ = bgc.PAR_half_saturation
    λ = bgc.PAR_attenuation_scale
    γ = bgc.dissolved_organic_phosphate_remin_timescale
    r = bgc.particulate_organic_phosphate_remin_timescale
    α = bgc.fraction_of_particulate_export
    Rᶜᴾ = bgc.stoichoimetric_ratio_carbon_to_phosphate   
    Rᴺᴾ = bgc.stoichoimetric_ratio_nitrate_to_phosphate  
    Rᴾᴼ = bgc.stoichoimetric_ratio_phosphate_to_oxygen   
    Rᶜᴺ = bgc.stoichoimetric_ratio_carbon_to_nitrate     
    Rᶜᴼ = bgc.stoichoimetric_ratio_carbon_to_oxygen      
    Rᶜᶠ = bgc.stoichoimetric_ratio_carbon_to_iron        
    Rᶠᴾ = bgc.stoichoimetric_ratio_iron_to_phosphate
    Lₜ     = bgc.ligand_concentration
    β     = bgc.ligand_stability_coefficient
    kˢᶜᵃᵛ = bgc.iron_scavenging_rate

    # Available photosynthetic radiation
    z = znode(i, j, k, grid, c, c, c)
    # TODO: design a user interface for prescribing incoming shortwave
    I = 700 * exp(z / λ)

    PO₄ = @inbounds fields.PO₄[i, j, k]
    NO₃ = @inbounds fields.NO₃[i, j, k]
    Feₜ = @inbounds fields.Fe[i, j, k]
    DOP = @inbounds fields.DOP[i, j, k]
    POP = @inbounds fields.POP[i, j, k]

    return Rᶠᴾ * (
                -   net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, PO₄, NO₃, Feₜ) 
                +   dissolved_organic_phosphate_remin(γ, DOP) 
                +   particulate_organic_phosphate_remin(r, POP)) 
                - iron_sources() 
                - iron_scavenging(kˢᶜᵃᵛ, Feₜ, Lₜ, β)
    end

"""
Tracer sources and sinks for DOP
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:DOP}, clock, fields)
    μᵖ = bgc.maximum_net_community_production_rate
    kᴺ = bgc.phosphate_half_saturation
    kᴾ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    kᴵ = bgc.PAR_half_saturation
    λ = bgc.PAR_attenuation_scale
    γ = bgc.dissolved_organic_phosphate_remin_timescale
    α = bgc.fraction_of_particulate_export     

    # Available photosynthetic radiation
    z = znode(i, j, k, grid, c, c, c)
    # TODO: design a user interface for prescribing incoming shortwave
    I = 700 * exp(z / λ)

    PO₄ = @inbounds fields.PO₄[i, j, k]
    NO₃ = @inbounds fields.NO₃[i, j, k]
    Feₜ = @inbounds fields.Fe[i, j, k]
    DOP = @inbounds fields.DOP[i, j, k]

    return - dissolved_organic_phosphate_remin(γ, DOP) +
             (1 - α) * net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, PO₄, NO₃, Feₜ)
end

"""
Tracer sources and sinks for POP
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:POP}, clock, fields)
    μᵖ = bgc.maximum_net_community_production_rate
    kᴺ = bgc.phosphate_half_saturation
    kᴾ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    kᴵ = bgc.PAR_half_saturation
    λ = bgc.PAR_attenuation_scale
    r = bgc.particulate_organic_phosphate_remin_timescale
    α = bgc.fraction_of_particulate_export     

    # Available photosynthetic radiation
    z = znode(i, j, k, grid, c, c, c)
    # TODO: design a user interface for prescribing incoming shortwave
    I = 700 * exp(z / λ)

    PO₄ = @inbounds fields.PO₄[i, j, k]
    NO₃ = @inbounds fields.NO₃[i, j, k]
    Feₜ = @inbounds fields.Fe[i, j, k]
    POP = @inbounds fields.POP[i, j, k]

    return α * net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, PO₄, NO₃, Feₜ) -
           particulate_organic_phosphate_remin(r, POP)
end