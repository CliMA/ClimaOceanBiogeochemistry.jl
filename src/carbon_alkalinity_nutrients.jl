using Oceananigans.Units: day
using Oceananigans.Grids: znode, Center
using Oceananigans.Fields: ZeroField, ConstantField
using Oceananigans.Biogeochemistry: AbstractBiogeochemistry

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers, biogeochemical_drift_velocity

const c = Center()

struct CarbonAlkalinityNutrients{FT} <: AbstractBiogeochemistry
    maximum_net_community_production_rate      :: FT # mol P m⁻³ s⁻¹
    phosphate_half_saturation                  :: FT # mol P m⁻³
    nitrate_half_saturation                    :: FT # mol N m⁻³
    iron_half_saturation                       :: FT # mol Fe m⁻³
    PAR_half_saturation                        :: FT  # W m⁻²
    PAR_attenuation_scale                      :: FT  # m
    fraction_of_particulate_export             :: FT
    dissolved_organic_phosphate_remin_timescale:: FT # s⁻¹
    stoichoimetric_ratio_carbon_to_phosphate   :: FT 
    stoichoimetric_ratio_nitrate_to_phosphate  :: FT 
    stoichoimetric_ratio_phosphate_to_oxygen   :: FT 
    stoichoimetric_ratio_iron_to_phosphate     :: FT 
    stoichoimetric_ratio_carbon_to_nitrate     :: FT 
    stoichoimetric_ratio_carbon_to_oxygen      :: FT 
    stoichoimetric_ratio_carbon_to_iron        :: FT 
    stoichoimetric_ratio_silicate_to_phosphate :: FT 
    rain_ratio_inorganic_to_organic_carbon     :: FT 
    martin_curve_exponent                      :: FT 
    iron_scavenging_rate                       :: FT # s⁻¹
    ligand_concentration                       :: FT # mol L m⁻³
    ligand_stability_coefficient               :: FT
end

"""
    CarbonAlkalinityNutrients(; reference_density                           = 1024,
                                maximum_net_community_production_rate       = 1 / day,
                                phosphate_half_saturation                   = 1e-7 * reference_density,
                                nitrate_half_saturation                     = 1.6e-6 * reference_density,
                                iron_half_saturation                        = 1e-10 * reference_density,
                                PAR_half_saturation                         = 10.0,
                                PAR_attenuation_scale                       = 25.0,
                                fraction_of_particulate_export              = 0.33
                                dissolved_organic_phosphate_remin_timescale = 1 / 30day,
                                stoichoimetric_ratio_carbon_to_phosphate    = 106.0
                                stoichoimetric_ratio_nitrate_to_phosphate   = 16.0
                                stoichoimetric_ratio_phosphate_to_oxygen    = 170.0,
                                stoichoimetric_ratio_phosphate_to_iron      = 4.68e-4
                                stoichoimetric_ratio_carbon_to_nitrate      = 106 / 16
                                stoichoimetric_ratio_carbon_to_oxygen       = 106 / 170,
                                stoichoimetric_ratio_carbon_to_iron         = 106 / 1.e-3
                                stoichoimetric_ratio_silicate_to_phosphate  = 15.0
                                rain_ratio_inorganic_to_organic_carbon      = 1e-1
                                martin_curve_exponent                       = 0.84,
                                iron_scavenging_rate                        = 5e-4 / day,
                                ligand_concentration                        = 1e-9 * reference_density,
                                ligand_stability_coefficient                = 1e8)

Return a six-tracer biogeochemistry model for the interaction of carbon, alkalinity, and nutrients.

Keyword Arguments
=================

Tracer names
============
* `DIC`: Dissolved Inorganic Carbon

* `Alk`: Alkalinity

* `PO₄`: Phosphate (macronutrient)

* `NO₃`: Nitrate (macronutrient)

* `DOP`: Dissolved Organic Phosphate (macronutrient)

* `Fe`: Dissolved Iron (micronutrient)

Biogeochemical functions
========================
* transitions for `DIC`, `Alk`, `PO₄`, `NO₃`, `DOP`, and `Fe`

* `biogeochemical_drift_velocity` for `D`, modeling the sinking of detritus at
  a constant `detritus_sinking_speed`.
"""
function CarbonAlkalinityNutrients(; reference_density = 1024,
                                   maximum_net_community_production_rate      = 1 / day, # mol P m⁻³ s⁻¹
                                   phosphate_half_saturation                  = 1e-7 * reference_density, # mol P m⁻³
                                   nitrate_half_saturation                    = 1.6e-6 * reference_density, # mol N m⁻³
                                   iron_half_saturation                       = 1e-10 * reference_density, # mol Fe m⁻³
                                   PAR_half_saturation                        = 10.0,  # W m⁻²
                                   PAR_attenuation_scale                      = 25.0,  # m
                                   fraction_of_particulate_export             = 0.33,
                                   dissolved_organic_phosphate_remin_timescale= 1 / 30day, # s⁻¹
                                   stoichoimetric_ratio_carbon_to_phosphate   = 106.0,
                                   stoichoimetric_ratio_nitrate_to_phosphate  = 16.0,
                                   stoichoimetric_ratio_phosphate_to_oxygen   = 170.0, 
                                   stoichoimetric_ratio_iron_to_phosphate     = 4.68e-4,
                                   stoichoimetric_ratio_carbon_to_nitrate     = 106 / 16,
                                   stoichoimetric_ratio_carbon_to_oxygen      = 106 / 170, 
                                   stoichoimetric_ratio_carbon_to_iron        = 106 / 1.e-3,
                                   stoichoimetric_ratio_silicate_to_phosphate = 15.0,
                                   rain_ratio_inorganic_to_organic_carbon     = 1e-1,
                                   martin_curve_exponent                      = 0.84,
                                   iron_scavenging_rate                       = 5e-4 / day, # s⁻¹
                                   ligand_concentration                       = 1e-9 * reference_density, # mol L m⁻³
                                   ligand_stability_coefficient               = 1e8)

    return CarbonAlkalinityNutrients(maximum_net_community_production_rate,
                                     phosphate_half_saturation,
                                     nitrate_half_saturation,
                                     iron_half_saturation,
                                     PAR_half_saturation,
                                     PAR_attenuation_scale,
                                     fraction_of_particulate_export,
                                     dissolved_organic_phosphate_remin_timescale,
                                     stoichoimetric_ratio_carbon_to_phosphate,
                                     stoichoimetric_ratio_nitrate_to_phosphate,
                                     stoichoimetric_ratio_phosphate_to_oxygen,
                                     stoichoimetric_ratio_iron_to_phosphate,
                                     stoichoimetric_ratio_carbon_to_nitrate,
                                     stoichoimetric_ratio_carbon_to_oxygen,
                                     stoichoimetric_ratio_carbon_to_iron,
                                     stoichoimetric_ratio_silicate_to_phosphate,
                                     rain_ratio_inorganic_to_organic_carbon,
                                     martin_curve_exponent,
                                     iron_scavenging_rate,
                                     ligand_concentration,
                                     ligand_stability_coefficient)
end

const CAN = CarbonAlkalinityNutrients

@inline required_biogeochemical_tracers(::CAN) = (:DIC, :Alk, :PO₄, :NO₃, :DOP, :Fe)

@inline function net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, P, N, F)
    return μᵖ * I / (I + kᴵ) * min(P / (P + kᴾ), N / (N + kᴺ), F / (F + kᶠ))
end

@inline dissolved_organic_phosphate_remin(γ, D) = γ * D

# Martin Curve
@inline particulate_organic_phosphate_remin() = 0.0

# exponential remineralization or below the lysocline
@inline particulate_inorganic_carbon_remin() = 0.0

@inline air_sea_flux_co2() = 0.0

@inline freshwater_virtual_flux() = 0.0

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

@inline iron_sources() = 0.0

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
    I = exp(z / λ)

    P = @inbounds fields.PO₄[i, j, k]
    N = @inbounds fields.NO₃[i, j, k]
    F = @inbounds fields.Fe[i, j, k]
    D = @inbounds fields.DOP[i, j, k]
    
    return Rᶜᴾ * (dissolved_organic_phosphate_remin(γ, D) -
                 (1 + Rᶜᵃᶜᵒ³ * α) * net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, P, N, F) +
                  particulate_organic_phosphate_remin()) +
           particulate_inorganic_carbon_remin() -
           air_sea_flux_co2() -
           freshwater_virtual_flux()
end

"""
Tracer sources and sinks for ALK
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:Alk}, clock, fields)
    μᵖ = bgc.maximum_net_community_production_rate
    kᴺ = bgc.phosphate_half_saturation
    kᴾ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    kᴵ = bgc.PAR_half_saturation
    λ = bgc.PAR_attenuation_scale
    γ = bgc.dissolved_organic_phosphate_remin_timescale
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
    I = exp(z / λ)

    P = @inbounds fields.PO₄[i, j, k]
    N = @inbounds fields.NO₃[i, j, k]
    F = @inbounds fields.Fe[i, j, k]
    D = @inbounds fields.DOP[i, j, k]
    
    return -Rᴺᴾ * (
        - (1 + Rᶜᵃᶜᵒ³ * α) * net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, P, N, F) +
        dissolved_organic_phosphate_remin(γ, D) +
        particulate_organic_phosphate_remin()) +
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
    I = exp(z / λ)

    P = @inbounds fields.PO₄[i, j, k]
    N = @inbounds fields.NO₃[i, j, k]
    F = @inbounds fields.Fe[i, j, k]
    D = @inbounds fields.DOP[i, j, k]

    return - net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, P, N, F) +
           dissolved_organic_phosphate_remin(γ, D) +
           particulate_organic_phosphate_remin()
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
    I = exp(z / λ)

    P = @inbounds fields.PO₄[i, j, k]
    N = @inbounds fields.NO₃[i, j, k]
    F = @inbounds fields.Fe[i, j, k]
    D = @inbounds fields.DOP[i, j, k]

    return Rᴺᴾ * (
           - net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, P, N, F) +
           dissolved_organic_phosphate_remin(γ, D) +
           particulate_organic_phosphate_remin())
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
    I = exp(z / λ)

    P = @inbounds fields.PO₄[i, j, k]
    N = @inbounds fields.NO₃[i, j, k]
    F = @inbounds fields.Fe[i, j, k]
    D = @inbounds fields.DOP[i, j, k]

    return (Rᶠᴾ * (
                - net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, P, N, F) +
                + dissolved_organic_phosphate_remin(γ, D)) 
           #    + particulate_organic_phosphate_remin()) +
           #    + iron_sources()
               - iron_scavenging(kˢᶜᵃᵛ, F, Lₜ, β))
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
    I = exp(z / λ)

    P = @inbounds fields.PO₄[i, j, k]
    N = @inbounds fields.NO₃[i, j, k]
    F = @inbounds fields.Fe[i, j, k]
    D = @inbounds fields.DOP[i, j, k]

    return - dissolved_organic_phosphate_remin(γ, D) +
             (1 - α) * net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, P, N, F)

end
