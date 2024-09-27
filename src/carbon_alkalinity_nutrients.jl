import Oceananigans.Biogeochemistry:
    biogeochemical_drift_velocity, required_biogeochemical_tracers

# using CUDA
using Adapt
import Adapt: adapt_structure, adapt
using Oceananigans.Biogeochemistry: AbstractBiogeochemistry
using Oceananigans.BoundaryConditions: ImpenetrableBoundaryCondition, fill_halo_regions!
using Oceananigans.Fields: ConstantField, ZeroField
using Oceananigans.Grids: Center, znode
using Oceananigans.Units: days

const c = Center()

struct CarbonAlkalinityNutrients{FT, W} <: AbstractBiogeochemistry
    reference_density                             :: FT
    maximum_net_community_production_rate         :: FT # mol PO₄ m⁻³ s⁻¹
    phosphate_half_saturation                     :: FT # mol PO₄ m⁻³
    nitrate_half_saturation                       :: FT # mol NO₃ m⁻³
    iron_half_saturation                          :: FT # mol Fe m⁻³
    incident_PAR                                  :: FT # W m⁻²
    PAR_half_saturation                           :: FT  # W m⁻²
    PAR_attenuation_scale                         :: FT  # m
    PAR_percent                                   :: FT  # m
    fraction_of_particulate_export                :: FT
    dissolved_organic_phosphorus_remin_timescale  :: FT # s⁻¹
    stoichoimetric_ratio_carbon_to_phosphate      :: FT 
    stoichoimetric_ratio_nitrate_to_phosphate     :: FT 
    stoichoimetric_ratio_phosphate_to_oxygen      :: FT 
    stoichoimetric_ratio_phosphate_to_iron        :: FT 
    stoichoimetric_ratio_carbon_to_nitrate        :: FT 
    stoichoimetric_ratio_carbon_to_oxygen         :: FT 
    stoichoimetric_ratio_carbon_to_iron           :: FT 
    stoichoimetric_ratio_silicate_to_phosphate    :: FT
    rain_ratio_inorganic_to_organic_carbon        :: FT 
    option_of_particulate_remin                   :: FT 
    particulate_organic_phosphate_remin_timescale :: FT
    iron_scavenging_rate                          :: FT # s⁻¹
    ligand_concentration                          :: FT # mol L m⁻³
    ligand_stability_coefficient                  :: FT
    martin_curve_exponent                         :: FT 
    particulate_organic_phosphorus_sinking_velocity   :: W  # m s⁻¹ 
end

"""
    CarbonAlkalinityNutrients(  reference_density                             = 1024.5,
                                maximum_net_community_production_rate         = 1 / day,
                                phosphate_half_saturation                     = 1e-7 * reference_density,
                                nitrate_half_saturation                       = 1.6e-6 * reference_density,
                                iron_half_saturation                          = 1e-10 * reference_density,
                                incident_PAR                                  = 700.0,
                                PAR_half_saturation                           = 10.0,
                                PAR_attenuation_scale                         = 25.0,
                                PAR_percent                                   = 0.01,
                                fraction_of_particulate_export                = 0.33,
                                dissolved_organic_phosphorus_remin_timescale   = 1 / 30day,
                                stoichoimetric_ratio_carbon_to_phosphate      = 106.0,
                                stoichoimetric_ratio_nitrate_to_phosphate     = 16.0,
                                stoichoimetric_ratio_phosphate_to_oxygen      = 170.0,
                                stoichoimetric_ratio_phosphate_to_iron        = 4.68e-4,
                                stoichoimetric_ratio_carbon_to_nitrate        = 106 / 16,
                                stoichoimetric_ratio_carbon_to_oxygen         = 106 / 170,
                                stoichoimetric_ratio_carbon_to_iron           = 106 / 1.e-3,
                                stoichoimetric_ratio_silicate_to_phosphate    = 15.0,
                                rain_ratio_inorganic_to_organic_carbon        = 1e-1,
                                option_of_particulate_remin                   = 1, 
                                iron_scavenging_rate                          = 5e-4 / day,
                                ligand_concentration                          = 1e-9 * reference_density,
                                ligand_stability_coefficient                  = 1e8,
                                martin_curve_exponent                         = 0.84,
                                particulate_organic_phosphorus_sinking_velocity   = -10.0 / day)

Return a seven-tracer biogeochemistry model for the interaction of carbon, alkalinity, and nutrients.

Keyword Arguments
=================

Tracer names
============
* `DIC`: Dissolved Inorganic Carbon

* `ALK`: Alkalinity

* `PO₄`: Phosphate (macronutrient)

* `NO₃`: Nitrate (macronutrient)

* `DOP`: Dissolved Organic Phosphorus (macronutrient)

* `POP`: Particulate Organic Phosphorus

* `Fe`: Dissolved Iron (micronutrient)

Biogeochemical functions
========================
* transitions for `DIC`, `ALK`, `PO₄`, `NO₃`, `DOP`, `POP` and `Fe`

* `biogeochemical_drift_velocity` for `POP`, modeling the sinking of POP at
  a constant `particulate_organic_phosphorus_sinking_velocity`.
"""
function CarbonAlkalinityNutrients(; grid,
                                   reference_density                            = 1024.5,
                                   maximum_net_community_production_rate        = 1 / day, # mol PO₄ m⁻³ s⁻¹
                                   phosphate_half_saturation                    = 1e-7 * reference_density, # mol PO₄ m⁻³
                                   nitrate_half_saturation                      = 1.6e-6 * reference_density, # mol NO₃ m⁻³
                                   iron_half_saturation                         = 1.e-10 * reference_density, # mol Fe m⁻³
                                   incident_PAR                                 = 700.0, # W m⁻²
                                   PAR_half_saturation                          = 10.0,  # W m⁻²
                                   PAR_attenuation_scale                        = 25.0,  # m
                                   PAR_percent                                  = 0.01,
                                   fraction_of_particulate_export               = 0.33,
                                   dissolved_organic_phosphorus_remin_timescale  = 1/30days, # s⁻¹
                                   stoichoimetric_ratio_carbon_to_phosphate     = 106.0,
                                   stoichoimetric_ratio_nitrate_to_phosphate    = 16.0,
                                   stoichoimetric_ratio_phosphate_to_oxygen     = 170.0, 
                                   stoichoimetric_ratio_phosphate_to_iron       = 4.68e-4,
                                   stoichoimetric_ratio_carbon_to_nitrate       = 106 / 16,
                                   stoichoimetric_ratio_carbon_to_oxygen        = 106 / 170, 
                                   stoichoimetric_ratio_carbon_to_iron          = 106 / 1e-3,
                                   stoichoimetric_ratio_silicate_to_phosphate   = 15.0,
                                   rain_ratio_inorganic_to_organic_carbon       = 1e-1,
                                   option_of_particulate_remin                  = 1, # r decrease with depth = 1; "power law" function = 2
                                   particulate_organic_phosphate_remin_timescale= 0.03 / day, 
                                   iron_scavenging_rate                         = 5e-4/days, # s⁻¹
                                   ligand_concentration                         = 1e-9 * reference_density, # mol L m⁻³
                                   ligand_stability_coefficient                 = 1e8,
                                   martin_curve_exponent                       = 0.84,
                                   particulate_organic_phosphorus_sinking_velocity  = -10.0 / day)

    if particulate_organic_phosphorus_sinking_velocity isa Number
            w₀ = particulate_organic_phosphorus_sinking_velocity
            no_penetration = ImpenetrableBoundaryCondition()
            bcs = FieldBoundaryConditions(grid, (Center, Center, Face),
                                        top=no_penetration, bottom=no_penetration)

            particulate_organic_phosphorus_sinking_velocity = ZFaceField(grid, boundary_conditions = bcs)

            set!(particulate_organic_phosphorus_sinking_velocity, w₀)

            fill_halo_regions!(particulate_organic_phosphorus_sinking_velocity)
    end
                            
    FT = eltype(grid)

    return CarbonAlkalinityNutrients(convert(FT, reference_density),
                                     convert(FT, maximum_net_community_production_rate),
                                     convert(FT, phosphate_half_saturation),
                                     convert(FT, nitrate_half_saturation),
                                     convert(FT, iron_half_saturation),
                                     convert(FT, incident_PAR),
                                     convert(FT, PAR_half_saturation),
                                     convert(FT, PAR_attenuation_scale),
                                     convert(FT, PAR_percent),
                                     convert(FT, fraction_of_particulate_export),
                                     convert(FT, dissolved_organic_phosphorus_remin_timescale),
                                     convert(FT, stoichoimetric_ratio_carbon_to_phosphate),
                                     convert(FT, stoichoimetric_ratio_nitrate_to_phosphate),
                                     convert(FT, stoichoimetric_ratio_phosphate_to_oxygen),
                                     convert(FT, stoichoimetric_ratio_phosphate_to_iron),
                                     convert(FT, stoichoimetric_ratio_carbon_to_nitrate),
                                     convert(FT, stoichoimetric_ratio_carbon_to_oxygen),
                                     convert(FT, stoichoimetric_ratio_carbon_to_iron),
                                     convert(FT, stoichoimetric_ratio_silicate_to_phosphate),
                                     convert(FT, rain_ratio_inorganic_to_organic_carbon),
                                     convert(FT, option_of_particulate_remin),
                                     convert(FT, particulate_organic_phosphate_remin_timescale), 
                                     convert(FT, iron_scavenging_rate),
                                     convert(FT, ligand_concentration),
                                     convert(FT, ligand_stability_coefficient),
                                     convert(FT, martin_curve_exponent),
                                     particulate_organic_phosphorus_sinking_velocity,
                                     )
end
    
const CAN = CarbonAlkalinityNutrients

Adapt.adapt_structure(to, bgc::CarbonAlkalinityNutrients) = 
    CarbonAlkalinityNutrients(adapt(to, bgc.reference_density),
                            adapt(to, bgc.maximum_net_community_production_rate),
                            adapt(to, bgc.phosphate_half_saturation),
                            adapt(to, bgc.nitrate_half_saturation),
                            adapt(to, bgc.iron_half_saturation),
                            adapt(to, bgc.incident_PAR),
                            adapt(to, bgc.PAR_half_saturation),
                            adapt(to, bgc.PAR_attenuation_scale),
                            adapt(to, bgc.PAR_percent),
                            adapt(to, bgc.fraction_of_particulate_export),
                            adapt(to, bgc.dissolved_organic_phosphorus_remin_timescale),
                            adapt(to, bgc.stoichoimetric_ratio_carbon_to_phosphate),
                            adapt(to, bgc.stoichoimetric_ratio_nitrate_to_phosphate),
                            adapt(to, bgc.stoichoimetric_ratio_phosphate_to_oxygen),
                            adapt(to, bgc.stoichoimetric_ratio_phosphate_to_iron),
                            adapt(to, bgc.stoichoimetric_ratio_carbon_to_nitrate),
                            adapt(to, bgc.stoichoimetric_ratio_carbon_to_oxygen),
                            adapt(to, bgc.stoichoimetric_ratio_carbon_to_iron),
                            adapt(to, bgc.stoichoimetric_ratio_silicate_to_phosphate),
                            adapt(to, bgc.rain_ratio_inorganic_to_organic_carbon),
                            adapt(to, bgc.option_of_particulate_remin),
                            adapt(to, bgc.particulate_organic_phosphate_remin_timescale), 
                            adapt(to, bgc.iron_scavenging_rate),
                            adapt(to, bgc.ligand_concentration),
                            adapt(to, bgc.ligand_stability_coefficient),
                            adapt(to, bgc.martin_curve_exponent),
                            adapt(to, bgc.particulate_organic_phosphorus_sinking_velocity))

@inline required_biogeochemical_tracers(::CAN) = (:DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe)

"""
Add a vertical sinking "drift velocity" for the particulate organic phosphorus (POP) tracer.
"""
@inline function biogeochemical_drift_velocity(bgc::CAN, ::Val{:POP})
    u = ZeroField()
    v = ZeroField()
    w = bgc.particulate_organic_phosphorus_sinking_velocity
    return (; u, v, w)
end

"""
    PAR(surface_photosynthetically_active_ratiation, 
        photosynthetically_active_ratiation_attenuation_scale, 
        depth)

Calculate the photosynthetically active radiation (PAR) at a given depth due to attenuation.
"""
@inline function PAR(surface_photosynthetically_active_ratiation, 
                     photosynthetically_active_ratiation_attenuation_scale, 
                     depth)

    I₀ = surface_photosynthetically_active_ratiation
    λ  = photosynthetically_active_ratiation_attenuation_scale
    z  = depth
    return I₀ * exp(z / λ)
end

"""
    net_community_production(maximum_net_community_production_rate,
                             light_half_saturation, 
                             phosphate_half_saturation, 
                             nitrate_half_saturation, 
                             iron_half_saturation, 
                             photosynthetically_active_radiation, 
                             phosphate_concentration, 
                             nitrate_concentration, 
                             iron_concentration)

"""
@inline function net_community_production(maximum_net_community_production_rate,
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
        P =phosphate_concentration
        N =nitrate_concentration
        F =iron_concentration

        # calculate the limitation terms
        Lₗᵢₘ = I / (I + kᴵ)
        Pₗᵢₘ = P / (P + kᴾ)
        Nₗᵢₘ = N / (N + kᴺ)
        Fₗᵢₘ = F / (F + kᶠ)

        # return the net community production
        return max(0,
        μᵖ * Lₗᵢₘ * min(Pₗᵢₘ, Nₗᵢₘ, Fₗᵢₘ)
        )
end

"""
dissolved_organic_phosphate_remin(remineralization_rate, 
                                  dissolved_organic_phosphorus_concentration)
Calculate the remineralization of dissolved organic phosphate.
"""
@inline dissolved_organic_phosphorus_remin(remineralization_rate, 
                                          dissolved_organic_phosphorus_concentration) = 
        max(0, remineralization_rate * dissolved_organic_phosphorus_concentration)

"""
Calculate remineralization of particulate organic phosphorus according to 
1) rate decreases with depth, 
or 2) a first-order rate constant .
"""
# @inline particulate_organic_phosphorus_remin(r, POP) = max(0, r * POP) # to avoid potential negative flux
@inline function particulate_organic_phosphorus_remin(option_of_particulate_remin,
                                                    particulate_organic_phosphate_remin_timescale, 
                                                    martin_curve_exponent,
                                                    particulate_organic_phosphorus_sinking_velocity,
                                                    PAR_attenuation_scale,
                                                    depth, percent_light,
                                                    particulate_organic_phosphorus_concentration)

        Rᵣ = option_of_particulate_remin                                  
        r = particulate_organic_phosphate_remin_timescale
        b = martin_curve_exponent    
        wₛ = particulate_organic_phosphorus_sinking_velocity
        λ = PAR_attenuation_scale
        z  = depth
        fᵢ= percent_light
        z₀ = log(fᵢ)*λ # The base of the euphotic layer depth (z₀) where PAR is degraded down to 1%     
        POP = particulate_organic_phosphorus_concentration

    return ifelse(Rᵣ == 1, max(0, b*wₛ/(z+z₀) * POP), max(0, r * POP))

end

# """
# Calculate the remineralization profile consistent with the Martin curve
# """
# struct ParticleReminParms{FT}
#     PAR_attenuation_scale   :: FT
#     martin_curve_exponent   :: FT
# end 

# @inline function martin_remin_profile(x, y, z, params::ParticleReminParms)
#     λ = params.PAR_attenuation_scale
#     b = params.martin_curve_exponent
#     z₀ = log(0.01)*λ 
#     return (min(z, z₀) / z₀)^-b
# end

"""
Calculate remineralization of particulate organic carbon.
"""
@inline particulate_inorganic_carbon_remin() = 0.0

@inline air_sea_flux_co2() = 0.0

@inline freshwater_virtual_flux() = 0.0

"""
    iron_scavenging(iron_scavenging_rate, iron_concentration, ligand_concentration, ligand_stability_coefficient)

Calculate the scavenging loss of iron. Iron scavenging depends on free iron, which involves solving a quadratic 
equation in terms of ligand concentration and stability coefficient. Ligand-complexed iron is safe from being scavenged.
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
       # Fₜ = FeL + Feᶠʳᵉᵉ
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
@inline iron_sources() = 0

"""
Tracer sources and sinks for Dissolved Inorganic Carbon (DIC)
"""

@inline function (bgc::CAN)(i, j, k, grid, ::Val{:DIC}, clock, fields)
    μᵖ = bgc.maximum_net_community_production_rate
    kᴺ = bgc.phosphate_half_saturation
    kᴾ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    I₀ = bgc.incident_PAR
    kᴵ = bgc.PAR_half_saturation
    λ = bgc.PAR_attenuation_scale
    fᵢ= bgc.PAR_percent
    γ = bgc.dissolved_organic_phosphorus_remin_timescale
    r = bgc.particulate_organic_phosphate_remin_timescale
    α = bgc.fraction_of_particulate_export
    Rᶜᴾ = bgc.stoichoimetric_ratio_carbon_to_phosphate       
    Rᶜᵃᶜᵒ³ = bgc.rain_ratio_inorganic_to_organic_carbon 
    b = bgc.martin_curve_exponent    
    wₛ = bgc.particulate_organic_phosphorus_sinking_velocity
    Rᵣ = bgc.option_of_particulate_remin

    # Available photosynthetic radiation
    z = znode(i, j, k, grid, c, c, c)
    I = PAR(I₀, λ, z)

    PO₄ = @inbounds fields.PO₄[i, j, k]
    NO₃ = @inbounds fields.NO₃[i, j, k]
    Feₜ = @inbounds fields.Fe[i, j, k]
    DOP = @inbounds fields.DOP[i, j, k]
    POP = @inbounds fields.POP[i, j, k]

    return (Rᶜᴾ * (
                    dissolved_organic_phosphorus_remin(γ, DOP) +
                    particulate_organic_phosphorus_remin(Rᵣ, r, b, wₛ[i,j,k], λ, z, fᵢ, POP) -
                    (1 + Rᶜᵃᶜᵒ³ * α) * net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, PO₄, NO₃, Feₜ)
                ) + particulate_inorganic_carbon_remin())
end

"""
Tracer sources and sinks for Alkalinity (ALK)
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:ALK}, clock, fields)
    μᵖ = bgc.maximum_net_community_production_rate
    kᴺ = bgc.phosphate_half_saturation
    kᴾ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    I₀ = bgc.incident_PAR
    kᴵ = bgc.PAR_half_saturation
    λ = bgc.PAR_attenuation_scale
    fᵢ= bgc.PAR_percent
    r = bgc.particulate_organic_phosphate_remin_timescale
    γ = bgc.dissolved_organic_phosphorus_remin_timescale
    α = bgc.fraction_of_particulate_export  
    Rᴺᴾ = bgc.stoichoimetric_ratio_nitrate_to_phosphate  
    Rᶜᵃᶜᵒ³ = bgc.rain_ratio_inorganic_to_organic_carbon     
    b = bgc.martin_curve_exponent    
    wₛ = bgc.particulate_organic_phosphorus_sinking_velocity
    Rᵣ = bgc.option_of_particulate_remin

    # Available photosynthetic radiation
    z = znode(i, j, k, grid, c, c, c)
    I = PAR(I₀, λ, z)

    PO₄ = @inbounds fields.PO₄[i, j, k]
    NO₃ = @inbounds fields.NO₃[i, j, k]
    Feₜ = @inbounds fields.Fe[i, j, k]
    DOP = @inbounds fields.DOP[i, j, k]
    POP = @inbounds fields.POP[i, j, k]

    return (-Rᴺᴾ * (
        - (1 + Rᶜᵃᶜᵒ³ * α) * net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, PO₄, NO₃, Feₜ) +
        dissolved_organic_phosphorus_remin(γ, DOP) +
        particulate_organic_phosphorus_remin(Rᵣ, r, b, wₛ[i,j,k], λ, z, fᵢ, POP)) +
        2 * particulate_inorganic_carbon_remin())
end

"""
Tracer sources and sinks for inorganic/dissolved Nitrate (NO₃).
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:NO₃}, clock, fields)
    μᵖ = bgc.maximum_net_community_production_rate
    kᴺ = bgc.phosphate_half_saturation
    kᴾ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    I₀ = bgc.incident_PAR
    kᴵ = bgc.PAR_half_saturation
    λ = bgc.PAR_attenuation_scale
    fᵢ= bgc.PAR_percent
    r = bgc.particulate_organic_phosphate_remin_timescale
    γ = bgc.dissolved_organic_phosphorus_remin_timescale 
    Rᴺᴾ = bgc.stoichoimetric_ratio_nitrate_to_phosphate  
    b = bgc.martin_curve_exponent    
    wₛ = bgc.particulate_organic_phosphorus_sinking_velocity
    α = bgc.fraction_of_particulate_export 
    Rᵣ = bgc.option_of_particulate_remin

    # Available photosynthetic radiation
    z = znode(i, j, k, grid, c, c, c)
    I = PAR(I₀, λ, z)

    PO₄ = @inbounds fields.PO₄[i, j, k]
    NO₃ = @inbounds fields.NO₃[i, j, k]
    Feₜ = @inbounds fields.Fe[i, j, k]
    DOP = @inbounds fields.DOP[i, j, k]
    POP = @inbounds fields.POP[i, j, k]

    return (Rᴺᴾ * (
           - net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, PO₄, NO₃, Feₜ) +
           dissolved_organic_phosphorus_remin(γ, DOP) +
           particulate_organic_phosphorus_remin(Rᵣ, r, b, wₛ[i,j,k], λ, z, fᵢ, POP)))
end

"""
Tracer sources and sinks for dissolved iron (FeT).
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:Fe}, clock, fields)
    μᵖ = bgc.maximum_net_community_production_rate
    kᴺ = bgc.phosphate_half_saturation
    kᴾ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    I₀ = bgc.incident_PAR
    kᴵ = bgc.PAR_half_saturation
    λ = bgc.PAR_attenuation_scale
    fᵢ= bgc.PAR_percent
    r = bgc.particulate_organic_phosphate_remin_timescale
    γ = bgc.dissolved_organic_phosphorus_remin_timescale      
    Rᶠᴾ = bgc.stoichoimetric_ratio_phosphate_to_iron
    Lₜ     = bgc.ligand_concentration
    β     = bgc.ligand_stability_coefficient
    kˢᶜᵃᵛ = bgc.iron_scavenging_rate
    b = bgc.martin_curve_exponent    
    wₛ = bgc.particulate_organic_phosphorus_sinking_velocity
    α = bgc.fraction_of_particulate_export 
    Rᵣ = bgc.option_of_particulate_remin

    # Available photosynthetic radiation
    z = znode(i, j, k, grid, c, c, c)
    I = PAR(I₀, λ, z)

    PO₄ = @inbounds fields.PO₄[i, j, k]
    NO₃ = @inbounds fields.NO₃[i, j, k]
    Feₜ = @inbounds fields.Fe[i, j, k]
    DOP = @inbounds fields.DOP[i, j, k]
    POP = @inbounds fields.POP[i, j, k]

    return (Rᶠᴾ * (
                -   net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, PO₄, NO₃, Feₜ) 
                +   dissolved_organic_phosphorus_remin(γ, DOP) 
                +   particulate_organic_phosphorus_remin(Rᵣ, r, b, wₛ[i,j,k], λ, z, fᵢ, POP)) +
            iron_sources() -
            iron_scavenging(kˢᶜᵃᵛ, Feₜ, Lₜ, β))
    end

    """
    Tracer sources and sinks for inorganic/dissolved phosphate (PO₄).
    """
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:PO₄}, clock, fields)
    μᵖ = bgc.maximum_net_community_production_rate
    kᴺ = bgc.phosphate_half_saturation
    kᴾ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    I₀ = bgc.incident_PAR
    kᴵ = bgc.PAR_half_saturation
    λ = bgc.PAR_attenuation_scale
    fᵢ= bgc.PAR_percent
    r = bgc.particulate_organic_phosphate_remin_timescale
    γ = bgc.dissolved_organic_phosphorus_remin_timescale 
    b = bgc.martin_curve_exponent    
    wₛ = bgc.particulate_organic_phosphorus_sinking_velocity
    α = bgc.fraction_of_particulate_export 
    Rᵣ = bgc.option_of_particulate_remin
    
    # Available photosynthetic radiation
    z = znode(i, j, k, grid, c, c, c)
    I = PAR(I₀, λ, z)
    
    PO₄ = @inbounds fields.PO₄[i, j, k]
    NO₃ = @inbounds fields.NO₃[i, j, k]
    Feₜ = @inbounds fields.Fe[i, j, k]
    DOP = @inbounds fields.DOP[i, j, k]
    POP = @inbounds fields.POP[i, j, k]

    return (- net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, PO₄, NO₃, Feₜ) +
            dissolved_organic_phosphorus_remin(γ, DOP) +
            particulate_organic_phosphorus_remin(Rᵣ, r, b, wₛ[i,j,k], λ, z, fᵢ, POP))
end

"""
Tracer sources and sinks for dissolved organic phosphorus (DOP).
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:DOP}, clock, fields)
    μᵖ = bgc.maximum_net_community_production_rate
    kᴺ = bgc.phosphate_half_saturation
    kᴾ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    I₀ = bgc.incident_PAR
    kᴵ = bgc.PAR_half_saturation
    λ = bgc.PAR_attenuation_scale
    γ = bgc.dissolved_organic_phosphorus_remin_timescale
    α = bgc.fraction_of_particulate_export     

    # Available photosynthetic radiation
    z = znode(i, j, k, grid, c, c, c)
    I = PAR(I₀, λ, z)

    PO₄ = @inbounds fields.PO₄[i, j, k]
    NO₃ = @inbounds fields.NO₃[i, j, k]
    Feₜ = @inbounds fields.Fe[i, j, k]
    DOP = @inbounds fields.DOP[i, j, k]

    return (- dissolved_organic_phosphorus_remin(γ, DOP) +
             (1 - α) * net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, PO₄, NO₃, Feₜ))
end

"""
Tracer sources and sinks for Particulate Organic Phosphorus (POP).
"""
@inline function (bgc::CAN)(i, j, k, grid, ::Val{:POP}, clock, fields)
    μᵖ = bgc.maximum_net_community_production_rate
    kᴺ = bgc.phosphate_half_saturation
    kᴾ = bgc.nitrate_half_saturation
    kᶠ = bgc.iron_half_saturation
    I₀ = bgc.incident_PAR
    kᴵ = bgc.PAR_half_saturation
    λ = bgc.PAR_attenuation_scale
    fᵢ= bgc.PAR_percent
    r = bgc.particulate_organic_phosphate_remin_timescale
    α = bgc.fraction_of_particulate_export     
    b = bgc.martin_curve_exponent    
    wₛ = bgc.particulate_organic_phosphorus_sinking_velocity
    Rᵣ = bgc.option_of_particulate_remin

    # Available photosynthetic radiation
    z = znode(i, j, k, grid, c, c, c)
    I = PAR(I₀, λ, z)

    PO₄ = @inbounds fields.PO₄[i, j, k]
    NO₃ = @inbounds fields.NO₃[i, j, k]
    Feₜ = @inbounds fields.Fe[i, j, k]
    POP = @inbounds fields.POP[i, j, k]

    return (α * net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, PO₄, NO₃, Feₜ) -
           particulate_organic_phosphorus_remin(Rᵣ, r, b, wₛ[i,j,k], λ, z, fᵢ, POP))
end