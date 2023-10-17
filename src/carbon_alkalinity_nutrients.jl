import Oceananigans.Biogeochemistry:
    biogeochemical_drift_velocity, required_biogeochemical_tracers

using Oceananigans.Biogeochemistry: AbstractBiogeochemistry
using Oceananigans.Fields: ConstantField, ZeroField, FunctionField
using Oceananigans.ImmersedBoundaries: active_cell
using Oceananigans.Grids: Center, Face, znode
using Oceananigans.Units: day
using Oceananigans.AbstractOperations: ∂z

"""
    CarbonAlkalinityNutrients(; kw...)

Return a six-tracer biogeochemistry model for the interaction of carbon, alkalinity, and nutrients.

Parameters
==========

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
struct CarbonAlkalinityNutrients{FT} <: AbstractBiogeochemistry
    maximum_net_community_production_rate       :: FT # mol P m⁻³ s⁻¹
    phosphate_half_saturation                   :: FT # mol P m⁻³
    nitrate_half_saturation                     :: FT # mol N m⁻³
    iron_half_saturation                        :: FT # mol Fe m⁻³
    PAR_half_saturation                         :: FT  # W m⁻²
    PAR_attenuation_scale                       :: FT  # m
    fraction_of_particulate_export              :: FT
    dissolved_organic_phosphate_remin_timescale :: FT # s⁻¹
    stoichoimetric_ratio_carbon_to_phosphate    :: FT 
    stoichoimetric_ratio_nitrate_to_phosphate   :: FT 
    stoichoimetric_ratio_phosphate_to_oxygen    :: FT 
    stoichoimetric_ratio_phosphate_to_iron      :: FT 
    stoichoimetric_ratio_carbon_to_nitrate      :: FT 
    stoichoimetric_ratio_carbon_to_oxygen       :: FT 
    stoichoimetric_ratio_carbon_to_iron         :: FT 
    stoichoimetric_ratio_silicate_to_phosphate  :: FT 
    rain_ratio_inorganic_to_organic_carbon      :: FT 
    martin_curve_exponent                       :: FT 
    iron_scavenging_rate                        :: FT # s⁻¹
    ligand_concentration                        :: FT # mol L m⁻³
    ligand_stability_coefficient                :: FT
    particulate_organic_phosphate_remin_exponent:: FT
    particulate_organic_phosphate_remin_curve   :: String
    particulate_organic_phosphate_bottom_remin  :: Bool
    particulate_inorganic_carbon_remin_exponent :: FT # s⁻¹
    particulate_inorganic_carbon_remin_curve    :: String
    particulate_inorganic_carbon_bottom_remin   :: Bool
end

function CarbonAlkalinityNutrients(; reference_density                         = 1024.5, # kg m⁻³
                                   maximum_net_community_production_rate       = 1 / day, # mol P m⁻³ s⁻¹
                                   phosphate_half_saturation                   = 1e-7 * reference_density, # mol P m⁻³
                                   nitrate_half_saturation                     = 1.6e-6 * reference_density, # mol N m⁻³
                                   iron_half_saturation                        = 1e-10 * reference_density, # mol Fe m⁻³
                                   PAR_half_saturation                         = 10.0,  # W m⁻²
                                   PAR_attenuation_scale                       = 25.0,  # m
                                   fraction_of_particulate_export              = 0.33,
                                   dissolved_organic_phosphate_remin_timescale = 1 / 30day, # s⁻¹
                                   stoichoimetric_ratio_carbon_to_phosphate    = 106.0, # mol C mol P⁻¹
                                   stoichoimetric_ratio_nitrate_to_phosphate   = 16.0, # mol N mol P⁻¹
                                   stoichoimetric_ratio_phosphate_to_oxygen    = 170.0, # mol P mol O⁻¹
                                   stoichoimetric_ratio_phosphate_to_iron      = 4.68e-4, # mol P mol Fe⁻¹
                                   stoichoimetric_ratio_carbon_to_nitrate      = 106 / 16, # mol C mol N⁻¹
                                   stoichoimetric_ratio_carbon_to_oxygen       = 106 / 170, # mol C mol O⁻¹
                                   stoichoimetric_ratio_carbon_to_iron         = 106 / 1.e-3, # mol C mol Fe⁻¹
                                   stoichoimetric_ratio_silicate_to_phosphate  = 15.0, # mol Si mol P⁻¹
                                   rain_ratio_inorganic_to_organic_carbon      = 1e-1, # mol C mol C⁻¹
                                   martin_curve_exponent                       = 0.84, #
                                   iron_scavenging_rate                        = 5e-4 / day, # s⁻¹
                                   ligand_concentration                        = 1e-9 * reference_density, # mol L m⁻³
                                   ligand_stability_coefficient                = 1e8, # s-1
                                   particulate_organic_phosphate_remin_exponent= 0.84,
                                   particulate_organic_phosphate_remin_curve   = "power law",
                                   particulate_organic_phosphate_bottom_remin  = false,
                                   particulate_inorganic_carbon_remin_exponent = 0.0,
                                   particulate_inorganic_carbon_remin_curve    = "exponential",
                                   particulate_inorganic_carbon_bottom_remin   = false,
                                   )

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
                                     stoichoimetric_ratio_phosphate_to_iron,
                                     stoichoimetric_ratio_carbon_to_nitrate,
                                     stoichoimetric_ratio_carbon_to_oxygen,
                                     stoichoimetric_ratio_carbon_to_iron,
                                     stoichoimetric_ratio_silicate_to_phosphate,
                                     rain_ratio_inorganic_to_organic_carbon,
                                     martin_curve_exponent,
                                     iron_scavenging_rate,
                                     ligand_concentration,
                                     ligand_stability_coefficient,
                                     particulate_organic_phosphate_remin_exponent,
                                     particulate_organic_phosphate_remin_curve,
                                     particulate_organic_phosphate_bottom_remin,
                                     particulate_inorganic_carbon_remin_exponent,
                                     particulate_inorganic_carbon_remin_curve,
                                     particulate_inorganic_carbon_bottom_remin,
                                     )
end

const CAN = CarbonAlkalinityNutrients

@inline required_biogeochemical_tracers(::CAN) = (:DIC, :Alk, :PO₄, :NO₃, :DOP, :Fe)

"""
    net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, P, N, F)
Calculate the rate of net community production, as the 
maximum growth rate, limited by the availability of light, phosphate, nitrate and iron.
"""
@inline function net_community_production(μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, P, N, F)
    return μᵖ .* (I ./ (I .+ kᴵ)) .* min.(P ./ (P .+ kᴾ), N ./ (N .+ kᴺ), F ./ (F .+ kᶠ))
end

"""
    dissolved_organic_phosphate_remin(γ, D)
Calculate the degradation of DOP to phosphate, with remineralization rate γ.
"""
@inline dissolved_organic_phosphate_remin(γ, D) = γ * D

struct ParticleReminParms{FT}
    remin_profile_reference_depth :: FT
    remin_profile_shape           :: String
    remin_profile_exponent        :: FT
end 

"""
    remin_profile(x, y, z, params::ParticleReminParms)
Calculate the remineralization profile, with parameters passed in a struct of type ParticleReminParms,
including remin_profile_reference_depth, remin_profile_shape, (either a "power law" or an "exponential"),
and any coefficients to the profile, such as remin_profile_exponent.
"""
@inline function remin_profile(x, y, z, params::ParticleReminParms)
    zʳᵉᶠ = params.remin_profile_reference_depth
    ɼ    = params.remin_profile_shape
    b    = params.remin_profile_exponent

    if ɼ == "power law"
        F = (min(z, zʳᵉᶠ) / zʳᵉᶠ)^-b
    elseif ɼ == "exponential"
        F = exp(-min(z, zʳᵉᶠ) / zʳᵉᶠ)
    end
    return F
end

"""
   particulate_remin(i, j, k, grid, fields, λ, α, μᵖ, kᴺ, kᴾ, kᶠ, kᴵ, ɼ = "power law", b = 0.84, ⎵ = false)
Calculate the degradation of sinking particles and return the tendency of local concentration. 
At depth z, there will be particles that  originate from the different surface layers, and 
remineralize according to a profile referenced to each layer. Thus we have to loop over all the 
productive surface layers to figure out how much remineralization affects local (deep) concentrations.

We can use light, depending on attenuation coefficient, λ, as a proxy otherwise would have to 
cycle throughout the entire watercolumn.

α, μᵖ, kᴺ, kᴾ, kᶠ, kᴵ relate to calculation of net community production in the productive surface layers. In
particular, α is the portion of production that is exported as sinking particles which could be set differently
for organic and inorganic particles such as organic carbon/phosphate or inorganic calcium carbonate.

ɼ is the remineralization profile, either "power law" or "exponential", b is the exponent of the remineralization 
profile (if it has one), and  and ⎵ ("the bottom") is a boolean that determines whether particles that reach 
the lowest active grid cell in depth are remineralized locally, or allowed  to export nutrients and carbon
out of the domain.
"""
@inline function particulate_remin(
    i, j, k, grid, 
    #α, μᵖ, kᴺ, kᴾ, kᶠ, kᴵ, Iʳᵉᶠ, Pʳᵉᶠ, Nʳᵉᶠ, Fʳᵉᶠ,
    PPʳᵉᶠ, 
    ɼ = "power law", b = 0.84, ⎵ = false
    )

    # Set the initial remineralization input to zero
    dPdt = 0.0

#    # Calculate the net community production at each reference depth above the current depth
#    PPʳᵉᶠ = α .* net_community_production(
#        μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, Iʳᵉᶠ, Pʳᵉᶠ, Nʳᵉᶠ, Fʳᵉᶠ
#        )

    # Loop over all layers above the current one, and determine if any particles
    #    are produced, and if so, how much remineralization there is.
    @inbounds for (ikʳᵉᶠ,kʳᵉᶠ) ∈ enumerate(k:grid.Nz)
    #if PPʳᵉᶠ[ikʳᵉᶠ] > zero(1.0) # or some small number?
    # There will be production of particles exported below this reference depth. 
    #   Calculate  depth change in the fraction of particle remineralization 
    #   given by the remin profile at current depth.
        zʳᵉᶠ = znode(i, j, kʳᵉᶠ, grid, Center(), Center(), Face())
        dFdz = ∂z(
            FunctionField{Center,Center,Face}(
                remin_profile,
                grid,
                clock      = nothing,
                parameters = ParticleReminParms(zʳᵉᶠ, ɼ, b),
            )
        )
    # ...then multiply by reference level export to calculate
    #  change in concentration due to remineralization of particles
        dPdt += dFdz[i, j, k] * PPʳᵉᶠ[ikʳᵉᶠ]
    
        # find if this is the bottom-most "active" grid cell and flux particles
        #  from this reference depth through the bottom of the domain?
        if ⎵ && ! active_cell(i,j,k-1,grid)
            # Integrate remin frac to find out what's left of the particles... 
            ∫F=Field(Integral(dFdz,dims=3))
            # ...and remineralize them here                
            dP_dz += (1 - ∫F[i, j]) * PPʳᵉᶠ[ikʳᵉᶠ]
        end
    #end
    end

    # Return total amount of inorganic constituent remineralized from sinking particles
    return dPdt
end

"""
    air_sea_flux_co2()
Calculate the air-sea flux of CO₂, according to the gas exchange parameterization 
of Wanninkhof (1992).
"""
@inline function air_sea_flux_co2()
    return 0.0
end

"""
    freshwater_virtual_flux()
Calculate the virtual flux of carbon and alkalinity, due to FW fluxes
"""
@inline function freshwater_virtual_flux()
    return 0.0
end

"""
    iron_scavenging(kˢᶜᵃᵛ, Fₜ, Lₜ, β)
Linear loss of free iron by scavenging onto sinking particles or precipitation.
"""
@inline function iron_scavenging(kˢᶜᵃᵛ, Fₜ, Lₜ, β)
       # solve for the equilibrium free iron concentration
       β⁻¹ = 1/β
       R₁  = 1
       R₂  = (Lₜ + β⁻¹ - Fₜ) 
       R₃  = -(Fₜ * β⁻¹) 

       # simple quadratic solution for roots
       discriminant = ( R₂*R₂ - ( 4*R₁*R₃ ))^(1/2)

       # directly solve for the free iron concentration
       Feᶠʳᵉᵉ = (-R₂ + discriminant) / (2*R₁) 

       # return the linear scavenging rate
       return - (kˢᶜᵃᵛ * Feᶠʳᵉᵉ)
end

@inline iron_sources() = 1e-6

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
    bᴾᴼᴾ = bgc.particulate_organic_phosphate_remin_exponent
    ɼᴾᴼᴾ = bgc.particulate_organic_phosphate_remin_curve
    ⎵ᴾᴼᴾ = bgc.particulate_organic_phosphate_bottom_remin
    bᴾᴵᶜ = bgc.particulate_inorganic_carbon_remin_exponent
    ɼᴾᴵᶜ = bgc.particulate_inorganic_carbon_remin_curve
    ⎵ᴾᴵᶜ = bgc.particulate_inorganic_carbon_bottom_remin

    # Available photosynthetic radiation
    z = znode(i, j, k:grid.Nz, grid, Center(), Center(), Center())
    
    I = exp.(z[1] ./ λ)
    P = @inbounds fields.PO₄[i, j, k]
    N = @inbounds fields.NO₃[i, j, k]
    F = @inbounds fields.Fe[ i, j, k]
    D = @inbounds fields.DOP[i, j, k]

    Iʳᵉᶠ = exp.(z ./ λ)
    Pʳᵉᶠ = @inbounds fields.PO₄[i, j, k:grid.Nz]
    Nʳᵉᶠ = @inbounds fields.NO₃[i, j, k:grid.Nz]
    Fʳᵉᶠ = @inbounds fields.Fe[ i, j, k:grid.Nz]
    Dʳᵉᶠ = @inbounds fields.DOP[i, j, k:grid.Nz]

    return Rᶜᴾ * (dissolved_organic_phosphate_remin(γ, D) -
                 (1 + Rᶜᵃᶜᵒ³ * α) * net_community_production(
                    μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, P, N, F
                    ) +
                 particulate_remin( 
                    i, j, k, grid, 
                    α.*net_community_production(
                        μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, Iʳᵉᶠ, Pʳᵉᶠ, Nʳᵉᶠ, Fʳᵉᶠ
                        ), 
                    ɼᴾᴼᴾ, bᴾᴼᴾ, ⎵ᴾᴼᴾ
                    )) +
                 #particulate_remin( i, j, k, grid, α,            μᵖ, kᴺ, kᴾ, kᶠ, kᴵ, Iʳᵉᶠ, Pʳᵉᶠ, Nʳᵉᶠ, Fʳᵉᶠ, ɼᴾᴼᴾ, bᴾᴼᴾ, ⎵ᴾᴼᴾ)) + #POC
                 #particulate_remin( i, j, k, grid, Rᶜᴾ*Rᶜᵃᶜᵒ³*α, μᵖ, kᴺ, kᴾ, kᶠ, kᴵ, Iʳᵉᶠ, Pʳᵉᶠ, Nʳᵉᶠ, Fʳᵉᶠ, ɼᴾᴵᶜ, bᴾᴵᶜ, ⎵ᴾᴵᶜ) -  #PIC
                 particulate_remin( 
                    i, j, k, grid, 
                    Rᶜᴾ.*Rᶜᵃᶜᵒ³.*α.*net_community_production(
                        μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, Iʳᵉᶠ, Pʳᵉᶠ, Nʳᵉᶠ, Fʳᵉᶠ
                        ), 
                    ɼᴾᴼᴾ, bᴾᴼᴾ, ⎵ᴾᴼᴾ
                    ) -
           - air_sea_flux_co2() -
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
    bᴾᴼᴾ = bgc.particulate_organic_phosphate_remin_exponent
    ɼᴾᴼᴾ = bgc.particulate_organic_phosphate_remin_curve
    ⎵ᴾᴼᴾ = bgc.particulate_organic_phosphate_bottom_remin
    bᴾᴵᶜ = bgc.particulate_inorganic_carbon_remin_exponent
    ɼᴾᴵᶜ = bgc.particulate_inorganic_carbon_remin_curve
    ⎵ᴾᴵᶜ = bgc.particulate_inorganic_carbon_bottom_remin

    # Available photosynthetic radiation
    z = znode(i, j, k:grid.Nz, grid, Center(), Center(), Center())

    I = exp.(z[1] ./ λ)
    P = @inbounds fields.PO₄[i, j, k]
    N = @inbounds fields.NO₃[i, j, k]
    F = @inbounds fields.Fe[ i, j, k]
    D = @inbounds fields.DOP[i, j, k]

    Iʳᵉᶠ = exp.(z ./ λ)
    Pʳᵉᶠ = @inbounds fields.PO₄[i, j, k:grid.Nz]
    Nʳᵉᶠ = @inbounds fields.NO₃[i, j, k:grid.Nz]
    Fʳᵉᶠ = @inbounds fields.Fe[ i, j, k:grid.Nz]
    Dʳᵉᶠ = @inbounds fields.DOP[i, j, k:grid.Nz]
    
    return -Rᴺᴾ * (
        - (1 + Rᶜᵃᶜᵒ³ * α) * net_community_production(
            μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, P, N, F
            ) +
        particulate_remin( 
            i, j, k, grid, 
            α.*net_community_production(
                μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, Iʳᵉᶠ, Pʳᵉᶠ, Nʳᵉᶠ, Fʳᵉᶠ
                ), 
            ɼᴾᴼᴾ, bᴾᴼᴾ, ⎵ᴾᴼᴾ
            ) +
    #    particulate_remin( i, j, k, grid, α,            μᵖ, kᴺ, kᴾ, kᶠ, kᴵ, Iʳᵉᶠ, Pʳᵉᶠ, Nʳᵉᶠ, Fʳᵉᶠ, ɼᴾᴼᴾ, bᴾᴼᴾ, ⎵ᴾᴼᴾ) + #POC
        dissolved_organic_phosphate_remin(γ, D)
        ) + 
     2 * particulate_remin(
        i, j, k, grid, 
        Rᶜᴾ.*Rᶜᵃᶜᵒ³.*α.*net_community_production(
            μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, Iʳᵉᶠ, Pʳᵉᶠ, Nʳᵉᶠ, Fʳᵉᶠ
            ), 
        ɼᴾᴼᴾ, bᴾᴼᴾ, ⎵ᴾᴼᴾ
        )   
    #2 * particulate_remin( i, j, k, grid, Rᶜᴾ*Rᶜᵃᶜᵒ³*α, μᵖ, kᴺ, kᴾ, kᶠ, kᴵ, Iʳᵉᶠ, Pʳᵉᶠ, Nʳᵉᶠ, Fʳᵉᶠ, ɼᴾᴵᶜ, bᴾᴵᶜ, ⎵ᴾᴵᶜ)  #PIC
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
    bᴾᴼᴾ = bgc.particulate_organic_phosphate_remin_exponent
    ɼᴾᴼᴾ = bgc.particulate_organic_phosphate_remin_curve
    ⎵ᴾᴼᴾ = bgc.particulate_organic_phosphate_bottom_remin

    # Available photosynthetic radiation
    z = znode(i, j, k:grid.Nz, grid, Center(), Center(), Center())
 
    I = exp.(z[1] ./ λ)
    P = @inbounds fields.PO₄[i, j, k]
    N = @inbounds fields.NO₃[i, j, k]
    F = @inbounds fields.Fe[ i, j, k]
    D = @inbounds fields.DOP[i, j, k]

    Iʳᵉᶠ = exp.(z ./ λ)
    Pʳᵉᶠ = @inbounds fields.PO₄[i, j, k:grid.Nz]
    Nʳᵉᶠ = @inbounds fields.NO₃[i, j, k:grid.Nz]
    Fʳᵉᶠ = @inbounds fields.Fe[ i, j, k:grid.Nz]

    return - net_community_production(
        μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, P, N, F
        ) +
#    particulate_remin( i, j, k, grid, α, μᵖ, kᴺ, kᴾ, kᶠ, kᴵ, Iʳᵉᶠ, Pʳᵉᶠ, Nʳᵉᶠ, Fʳᵉᶠ, ɼᴾᴼᴾ, bᴾᴼᴾ, ⎵ᴾᴼᴾ) +
    particulate_remin(
        i, j, k, grid, 
        α.*net_community_production(
            μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, Iʳᵉᶠ, Pʳᵉᶠ, Nʳᵉᶠ, Fʳᵉᶠ
            ), 
        ɼᴾᴼᴾ, bᴾᴼᴾ, ⎵ᴾᴼᴾ
        ) +
    dissolved_organic_phosphate_remin(γ, D)
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
    bᴾᴼᴾ = bgc.particulate_organic_phosphate_remin_exponent
    ɼᴾᴼᴾ = bgc.particulate_organic_phosphate_remin_curve
    ⎵ᴾᴼᴾ = bgc.particulate_organic_phosphate_bottom_remin

    # Available photosynthetic radiation
    z = znode(i, j, k:grid.Nz, grid, Center(), Center(), Center())

    I = exp.(z[1] ./ λ)
    P = @inbounds fields.PO₄[i, j, k]
    N = @inbounds fields.NO₃[i, j, k]
    F = @inbounds fields.Fe[ i, j, k]
    D = @inbounds fields.DOP[i, j, k]

    Iʳᵉᶠ = exp.(z ./ λ)
    Pʳᵉᶠ = @inbounds fields.PO₄[i, j, k:grid.Nz]
    Nʳᵉᶠ = @inbounds fields.NO₃[i, j, k:grid.Nz]
    Fʳᵉᶠ = @inbounds fields.Fe[ i, j, k:grid.Nz]
    Dʳᵉᶠ = @inbounds fields.DOP[i, j, k:grid.Nz]

    return Rᴺᴾ * (
           - net_community_production(
            μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, P, N, F
            ) +
           particulate_remin(
                i, j, k, grid, α.*net_community_production(
                    μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, Iʳᵉᶠ, Pʳᵉᶠ, Nʳᵉᶠ, Fʳᵉᶠ
                    ), 
                ɼᴾᴼᴾ, bᴾᴼᴾ, ⎵ᴾᴼᴾ
                ) +
           #particulate_remin( i, j, k, grid, α, μᵖ, kᴺ, kᴾ, kᶠ, kᴵ, Iʳᵉᶠ, Pʳᵉᶠ, Nʳᵉᶠ, Fʳᵉᶠ, ɼᴾᴼᴾ, bᴾᴼᴾ, ⎵ᴾᴼᴾ) + #POP
           dissolved_organic_phosphate_remin(γ, D)
           )
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
    Rᶠᴾ = bgc.stoichoimetric_ratio_phosphate_to_iron
    Lₜ = bgc.ligand_concentration
    β = bgc.ligand_stability_coefficient
    kˢᶜᵃᵛ = bgc.iron_scavenging_rate
    bᴾᴼᴾ = bgc.particulate_organic_phosphate_remin_exponent
    ɼᴾᴼᴾ = bgc.particulate_organic_phosphate_remin_curve
    ⎵ᴾᴼᴾ = bgc.particulate_organic_phosphate_bottom_remin

    # Available photosynthetic radiation
    z = znode(i, j, k:grid.Nz, grid, Center(), Center(), Center())

    I = exp.(z[1] ./ λ)
    P = @inbounds fields.PO₄[i, j, k]
    N = @inbounds fields.NO₃[i, j, k]
    F = @inbounds fields.Fe[ i, j, k]
    D = @inbounds fields.DOP[i, j, k]

    Iʳᵉᶠ = exp.(z ./ λ)
    Pʳᵉᶠ = @inbounds fields.PO₄[i, j, k:grid.Nz]
    Nʳᵉᶠ = @inbounds fields.NO₃[i, j, k:grid.Nz]
    Fʳᵉᶠ = @inbounds fields.Fe[ i, j, k:grid.Nz]
    Dʳᵉᶠ = @inbounds fields.DOP[i, j, k:grid.Nz]

    return Rᶠᴾ * (
                - net_community_production(
                    μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, P, N, F
                    ) +
                 dissolved_organic_phosphate_remin(γ, D) +
                 particulate_remin( 
                    i, j, k, grid, α.*net_community_production(
                        μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, Iʳᵉᶠ, Pʳᵉᶠ, Nʳᵉᶠ, Fʳᵉᶠ
                        ), 
                    ɼᴾᴼᴾ, bᴾᴼᴾ, ⎵ᴾᴼᴾ
                    ) +
                 #particulate_remin( i, j, k, grid, α, μᵖ, kᴺ, kᴾ, kᶠ, kᴵ, Iʳᵉᶠ, Pʳᵉᶠ, Nʳᵉᶠ, Fʳᵉᶠ, ɼᴾᴼᴾ, bᴾᴼᴾ, ⎵ᴾᴼᴾ) + #POC
                 iron_sources() +
                 iron_scavenging(kˢᶜᵃᵛ, F, Lₜ, β))
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
    z = znode(i, j, k:grid.Nz, grid, Center(), Center(), Center())

    I = exp.(z[1] ./ λ)
    P = @inbounds fields.PO₄[i, j, k]
    N = @inbounds fields.NO₃[i, j, k]
    F = @inbounds fields.Fe[ i, j, k]
    D = @inbounds fields.DOP[i, j, k]

    Iʳᵉᶠ = exp.(z ./ λ)
    Pʳᵉᶠ = @inbounds fields.PO₄[i, j, k:grid.Nz]
    Nʳᵉᶠ = @inbounds fields.NO₃[i, j, k:grid.Nz]
    Fʳᵉᶠ = @inbounds fields.Fe[ i, j, k:grid.Nz]
    Dʳᵉᶠ = @inbounds fields.DOP[i, j, k:grid.Nz]

    return - dissolved_organic_phosphate_remin(γ, D) +
             (1 - α) * net_community_production(
                μᵖ, kᴵ, kᴾ, kᴺ, kᶠ, I, P, N, F
                )
end
