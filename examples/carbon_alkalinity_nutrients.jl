using Oceananigans.Units: day
using Oceananigans.Grids: znode, Center
using Oceananigans.Fields: ZeroField, ConstantField
using Oceananigans.Biogeochemistry: AbstractBiogeochemistry

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers, biogeochemical_drift_velocity

const c = Center()

"""
    CarbonAlkalinityNutrients(; kw...)

Return a five-tracer biogeochemistry model for the interaction of carbon, alkalinity, and nutrients.

Parameters
==========
  * `maximum_net_community_production_rate`: (mol P m⁻³ s⁻¹) Net community production rate unlimited by the
                                     availability of nutrients and light. Default: 1/day.

  * `maximum_bacteria_growth_rate`: (s⁻¹) Growth rate of plankton `B` unlimited by the
                                    availability of nutrients and light. Default = 0.5/day.

  * `bacteria_yield`: Determines fractional nutrient production by bacteria production
                      relative to consumption of detritus such that ``∂_t N / ∂_t D = 1 - y``,
                      where `y = bacteria_yield`. Default: 0.2.

  * `quadratic_mortality_rate`: (s⁻¹) Mortality rate of both plankton and bacteria.

  * `nutrient_half_saturation`: (mmol m⁻³) Half-saturation of nutrients for plankton production.

  * `detritus_half_saturation`: (mmol m⁻³) Half-saturation of nutrients for bacteria production.
                                Deafult = 10.0 mmol m⁻³.

  * `PAR_half_saturation`: (W m⁻²) Half-saturation of photosynthetically available radiation (PAR)
                           for plankton production.

  * `PAR_attenuation_scale`: (m) Depth scale over which photosynthetically available radiation (PAR)
                             attenuates exponentially.

  * `detritus_sinking_speed`: (m s⁻¹) Sinking velocity of detritus.

Tracer names
============
  * `C`: Dissolved Inorganic Carbon
  * `A`: Alkalinity
  * `P`: Phosphate
  * `N`: Nitrate 
  * `D`: Dissolved Organic Phosphate
  * `F`: Dissolved Iron

Biogeochemical functions
========================
  * transitions for `C`, `A`, `P`, `N`, `D`, and `F`
  * `biogeochemical_drift_velocity` for `D`, modeling the sinking of detritus at
    a constant `detritus_sinking_speed`.
"""
const conv_molkg_molm3 = 1024.5
Base.@kwdef struct CarbonAlkalinityNutrients{FT} <: AbstractBiogeochemistry
    maximum_net_community_production_rate      :: FT = 1/day # mol P m⁻³ s⁻¹
    phosphate_half_saturation                  :: FT = 0.1e-6*conv_molkg_molm3 # mol P m⁻³
    nitrate_half_saturation                    :: FT = 16.*0.1e-6*conv_molkg_molm3 # mol N m⁻³
    iron_half_saturation                       :: FT = 0.1e-9*conv_molkg_molm3 # mol Fe m⁻³
    PAR_half_saturation                        :: FT = 10.0  # W m⁻²
    PAR_attenuation_scale                      :: FT = 25.0  # m
    fraction_of_particulate_export             :: FT = 0.33
    dissolved_organic_phosphate_remin_timescale:: FT = 1/(30*day) # s⁻¹
    stoichoimetric_ratio_carbon_to_phosphate   :: FT = 106.
    stoichoimetric_ratio_nitrate_to_phosphate  :: FT = 16.
    stoichoimetric_ratio_phosphate_to_oxygen   :: FT = 170. 
    stoichoimetric_ratio_phosphate_to_iron     :: FT = 4.68e-4    
    stoichoimetric_ratio_carbon_to_nitrate     :: FT = 106./16. 
    stoichoimetric_ratio_carbon_to_oxygen      :: FT = 106./170. 
    stoichoimetric_ratio_carbon_to_iron        :: FT = 106./1.e-3 
    stoichoimetric_ratio_silicate_to_phosphate :: FT = 15.
    rain_ratio_inorganic_to_organic_carbon     :: FT = 10.e-2 # %
    martin_curve_exponent                      :: FT = 0.84
    iron_scavenging_rate                       :: FT = 5.e-4/day # s⁻¹
    ligand_concentration                       :: FT = 1.e-9*conv_molkg_molm3 # mol L m⁻³
    ligand_stability_coefficient               :: FT = 1.e8
end

const CAN = CarbonAlkalinityNutrients

@inline required_biogeochemical_tracers(::CAN) = (:C, :A, :P, :N, :D, :F)

# For implicit time stepping, consider:
#
# ∂t P ≈ (P⁺ - P⁻) / Δt = - m * P⁻
#   -> P⁺ = P⁻ - Δt * m * P⁻
#
# (P⁺ - P⁻) / Δt = - m * P⁺
#   -> P⁺ = P⁻ / (1 + Δt * m)
#

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
@inline iron_scavenging(kˢᶜᵃᵛ, F, Lₜ, β) = kˢᶜᵃᵛ * ( F - Lₜ )

@inline iron_sources() = 0.0
"""
Tracer sources and sinks for DIC
"""
@inline function (bgc::NutrientsPlanktonBacteriaDetritus)(i, j, k, grid, ::Val{:C}, clock, fields)
    μᵖ = bgc.maximum_plankton_growth_rate
    kᴰ = bgc.detritus_half_saturation
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

    P = @inbounds fields.P[i, j, k]
    N = @inbounds fields.N[i, j, k]
    F = @inbounds fields.F[i, j, k]
    D = @inbounds fields.D[i, j, k]
    
    return Rᶜᴾ * (
           - (1 + Rᶜᵃᶜᵒ³ * α) * net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, P, N, F) \
           + dissolved_organic_phosphate_remin(γ, D) \
           + particulate_organic_phosphate_remin() \
           ) \
           + particulate_inorganic_carbon_remin() \
           - air_sea_flux_co2()
           - freshwater_virtual_flux()
end

"""
Tracer sources and sinks for ALK
"""
@inline function (bgc::NutrientsPlanktonBacteriaDetritus)(i, j, k, grid, ::Val{:C}, clock, fields)
    μᵖ = bgc.maximum_plankton_growth_rate
    kᴰ = bgc.detritus_half_saturation
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

    P = @inbounds fields.P[i, j, k]
    N = @inbounds fields.N[i, j, k]
    F = @inbounds fields.F[i, j, k]
    D = @inbounds fields.D[i, j, k]
    
    return -Rᴺᴾ * (
        - (1 + Rᶜᵃᶜᵒ³ * α) * net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, P, N, F) \
        + dissolved_organic_phosphate_remin(γ, D) \
        + particulate_organic_phosphate_remin() \
        ) \
        + 2.0 * particulate_inorganic_carbon_remin()
end

"""
Tracer sources and sinks for PO₄
"""
@inline function (bgc::NutrientsPlanktonBacteriaDetritus)(i, j, k, grid, ::Val{:P}, clock, fields)
    μᵖ = bgc.maximum_plankton_growth_rate
    kᴰ = bgc.detritus_half_saturation
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

    P = @inbounds fields.P[i, j, k]
    N = @inbounds fields.N[i, j, k]
    F = @inbounds fields.F[i, j, k]
    D = @inbounds fields.D[i, j, k]

    return - net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, P, N, F) \
           + dissolved_organic_phosphate_remin(γ, D) \
           + particulate_organic_phosphate_remin()
end

"""
Tracer sources and sinks for NO₃
"""
@inline function (bgc::NutrientsPlanktonBacteriaDetritus)(i, j, k, grid, ::Val{:N}, clock, fields)
    μᵖ = bgc.maximum_plankton_growth_rate
    kᴰ = bgc.detritus_half_saturation
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

    P = @inbounds fields.P[i, j, k]
    N = @inbounds fields.N[i, j, k]
    F = @inbounds fields.F[i, j, k]
    D = @inbounds fields.D[i, j, k]

    return Rᴺᴾ * ( \
           - net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, P, N, F) \
           + dissolved_organic_phosphate_remin(γ, D) \
           + particulate_organic_phosphate_remin() \
           )
end

"""
Tracer sources and sinks for FeT
"""
@inline function (bgc::NutrientsPlanktonBacteriaDetritus)(i, j, k, grid, ::Val{:F}, clock, fields)
    μᵖ = bgc.maximum_plankton_growth_rate
    kᴰ = bgc.detritus_half_saturation
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

    P = @inbounds fields.P[i, j, k]
    N = @inbounds fields.N[i, j, k]
    F = @inbounds fields.F[i, j, k]
    D = @inbounds fields.D[i, j, k]

    return Rᶠᴾ * ( \
           - net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, P, N, F) \
           + dissolved_organic_phosphate_remin(γ, D) \
           + particulate_organic_phosphate_remin() \
           ) \
           + iron_sources()
           - iron_scavenging()
    end

"""
Tracer sources and sinks for DOP
"""
@inline function (bgc::NutrientsPlanktonBacteriaDetritus)(i, j, k, grid, ::Val{:D}, clock, fields)
    μᵖ = bgc.maximum_plankton_growth_rate
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

    P = @inbounds fields.P[i, j, k]
    N = @inbounds fields.N[i, j, k]
    F = @inbounds fields.F[i, j, k]
    D = @inbounds fields.D[i, j, k]

    return - dissolved_organic_phosphate_remin(γ, D)
           + (1 - α) * net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, P, N, F)

end
