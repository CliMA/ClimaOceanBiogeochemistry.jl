using Oceananigans
using Oceananigans.Units: day
using Oceananigans.Grids: znode, Center, AbstractTopology, Flat, Bounded
using Oceananigans.BoundaryConditions: ImpenetrableBoundaryCondition, fill_halo_regions!
using Oceananigans.Fields: ZeroField, ZFaceField
using Oceananigans.Biogeochemistry: AbstractBiogeochemistry

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers, biogeochemical_drift_velocity

const c = Center()

struct NutrientsPlanktonBacteriaDetritus{FT, W} <: AbstractBiogeochemistry
    maximum_plankton_growth_rate :: FT 
    maximum_bacteria_growth_rate :: FT 
    maximum_grazing_rate :: FT         
    bacteria_yield :: FT               
    zooplankton_yield :: FT            
    linear_remineralization_rate :: FT 
    linear_mortality_rate :: FT        
    quadratic_mortality_rate :: FT     
    quadratic_mortality_rate_Z :: FT   
    nutrient_half_saturation :: FT     
    detritus_half_saturation :: FT     
    grazing_half_saturation  :: FT     
    PAR_half_saturation :: FT          
    PAR_attenuation_scale :: FT        
    detritus_vertical_velocity :: W        
end

"""
    NutrientsPlanktonBacteriaDetritus(; grid,
                                        maximum_plankton_growth_rate = 1/day,
                                        maximum_bacteria_growth_rate = 1/day
                                        maximum_grazing_rate         = 3/day
                                        bacteria_yield               = 0.2
                                        zooplankton_yield            = 0.3
                                        linear_remineralization_rate = 0.03/day,
                                        linear_mortality_rate        = 0.01/day,
                                        quadratic_mortality_rate     = 0.1/day,
                                        quadratic_mortality_rate_Z   = 1/day,
                                        nutrient_half_saturation     = 0.1,
                                        detritus_half_saturation     = 0.1,
                                        grazing_half_saturation      = 3.0,
                                        PAR_half_saturation          = 10.0,
                                        PAR_attenuation_scale        = 25.0,
                                        detritus_vertical_velocity   = -10/day)

Return a six-tracer biogeochemistry model for the interaction of nutrients (N), phytoplankton (P), 
zooplankton(Z), bacteria (B), dissolved detritus (D1), and particulate detritus (D2).

Keyword Arguments
=================
* `grid` (required): An Oceananigans' grid.

* `maximum_plankton_growth_rate`: (s⁻¹) Growth rate of plankton `P` unlimited by the
                                    availability of nutrients and light. Default: 1/day.

* `maximum_bacteria_growth_rate`: (s⁻¹) Growth rate of plankton `B` unlimited by the
                                  availability of nutrients and light. Default = 0.5/day.

* `maximum_grazing_rate`: (s⁻¹) Maximum grazing rate of phytoplankton by zooplankton.

* `bacteria_yield`: Determines fractional nutrient production by bacteria production 
                    relative to consumption of detritus such that ``∂_t N / ∂_t D = 1 - y``,
                    where `y = bacteria_yield`. Default: 0.2.

* `linear_remineralization_rate`: (s⁻¹) Remineralization rate constant of detritus 'D', 
                                  assuming linear remineralization of 'D', while 
                                  implicitly modeling bacteria 'B'. Default = 0.3/day.

* `linear_mortality_rate`: (s⁻¹) Linear term of the mortality rate of both plankton and bacteria.

* `quadratic_mortality_rate`: (s⁻¹) Quadratic term of the mortality rate of both plankton and bacteria.

* `nutrient_half_saturation`: (mmol m⁻³) Half-saturation of nutrients for plankton production.

* `detritus_half_saturation`: (mmol m⁻³) Half-saturation of nutrients for bacteria production.
                              Default = 10.0 mmol m⁻³.

* `phytoplankton_half_saturation`: (mmol m⁻³) Half-saturation of phytoplankton for zooplankton production.

* `zooplankton_assimilation`: Fractional assimilation efficiency for zooplankton.

* `PAR_half_saturation`: (W m⁻²) Half-saturation of photosynthetically available radiation (PAR)
                         for plankton production.

* `PAR_attenuation_scale`: (m) Depth scale over which photosynthetically available radiation (PAR)
                            attenuates exponentially.

* `detritus_sinking_speed`: (m s⁻¹) Sinking velocity of particulate detritus.

Tracer names
============
* `N`: nutrients

* `P`: phytoplankton

* `Z`: zooplankton

* `B`: bacteria

* `D1`: detritus 1 - dissolved

* `D2`: detritus 2 - particulate

Biogeochemical functions
========================
* transitions for `N`, `P`, `Z`, `B`, `D1`, `D2`

* `biogeochemical_drift_velocity` for `D2`, modeling the sinking of detritus at
  a constant `detritus_sinking_speed`.
"""
function NutrientsPlanktonBacteriaDetritus(; grid,
                                           maximum_plankton_growth_rate = 1/day, # Add reference for each parameter
                                           maximum_bacteria_growth_rate = 1/day,
                                           maximum_grazing_rate         = 3/day,
                                           bacteria_yield               = 0.2,
                                           zooplankton_yield            = 0.3,
                                           linear_remineralization_rate = 0.03/day, 
                                           linear_mortality_rate        = 0.01/day, # m³/mmol/day
                                           quadratic_mortality_rate     = 0.1/day,  # m³/mmol/day
                                           quadratic_mortality_rate_Z   = 1/day,    # m³/mmol/day (zooplankton quadratic mortality)
                                           nutrient_half_saturation     = 0.1,      # mmol m⁻³
                                           detritus_half_saturation     = 0.1,      # mmol m⁻³
                                           grazing_half_saturation      = 3.0,      # mmol m⁻³
                                           PAR_half_saturation          = 10.0,     # W m⁻²
                                           PAR_attenuation_scale        = 25.0,     # m
                                           detritus_vertical_velocity   = -10/day)  # m s⁻¹

    if detritus_vertical_velocity isa Number
        w₀ = detritus_vertical_velocity
        no_penetration = ImpenetrableBoundaryCondition()

        bcs = FieldBoundaryConditions(grid, (Center, Center, Face),
                                      top=no_penetration, bottom=no_penetration)

        detritus_vertical_velocity = ZFaceField(grid, boundary_conditions = bcs)

        set!(detritus_vertical_velocity, w₀)

        fill_halo_regions!(detritus_vertical_velocity)
    end

    FT = eltype(grid)

    return NutrientsPlanktonBacteriaDetritus(convert(FT, maximum_plankton_growth_rate),   
                                             convert(FT, maximum_bacteria_growth_rate),   
                                             convert(FT, maximum_grazing_rate),           
                                             convert(FT, bacteria_yield),                 
                                             convert(FT, zooplankton_yield),              
                                             convert(FT, linear_remineralization_rate),   
                                             convert(FT, linear_mortality_rate),          
                                             convert(FT, quadratic_mortality_rate),       
                                             convert(FT, quadratic_mortality_rate_Z),     
                                             convert(FT, nutrient_half_saturation),       
                                             convert(FT, detritus_half_saturation),       
                                             convert(FT, grazing_half_saturation),        
                                             convert(FT, PAR_half_saturation),            
                                             convert(FT, PAR_attenuation_scale),          
                                             detritus_vertical_velocity)
end

const NPZBD = NutrientsPlanktonBacteriaDetritus

@inline required_biogeochemical_tracers(::NPZBD) = (:N, :P, :Z, :B, :D1, :D2)

@inline function biogeochemical_drift_velocity(bgc::NPZBD, ::Val{:D2})
    u = ZeroField()
    v = ZeroField()
    w = bgc.detritus_vertical_velocity
    return (; u, v, w)
end

# For implicit time stepping, consider:
#
# ∂t P ≈ (P⁺ - P⁻) / Δt = - m * P⁻
#   -> P⁺ = P⁻ - Δt * m * P⁻
#
# (P⁺ - P⁻) / Δt = - m * P⁺
#   -> P⁺ = P⁻ / (1 + Δt * m)
#

# A depth-dependent temperature curve from Zakem (2018)
# Temp = 12 .*exp.(z./ 150) .+ 12 .*exp.(z ./ 500) .+ 2

# Temperature modification to metabolic rates, following the Arrhenius equation
# @inline temp_fun(Temp) = 0.8 .* exp.(-4000 .*(1 ./ (Temp .+ 273.15) .- 1 ./ 293.15))

@inline bacteria_production(μᵇ, kᴰ, y, D, B) = y * μᵇ * D / (D + kᴰ) * B 
@inline phytoplankton_production(μᵖ, kᴺ, kᴵ, I, N, P) = (μᵖ * min(N / (N + kᴺ) , I / (I + kᴵ)) * P) 
@inline zooplankton_graze_phytoplankton(gₘ, kᵍ, γ, P, Z) = γ * gₘ * P / (P + kᵍ) * Z
@inline zooplankton_graze_bacteria(gₘ, kᵍ, γ, B, Z) = γ * gₘ * B / (B + kᵍ) * Z
@inline bacteria_mortality(mlin,mq, B) = mlin * B + mq * B^2
@inline phytoplankton_mortality(mlin,mq, P) = mlin * P + mq * P^2
@inline zooplankton_mortality(mlin,mq_Z, Z) = mlin * Z + mq_Z * Z^2
@inline detritus_remineralization(r, D) = r * D

@inline function (bgc::NutrientsPlanktonBacteriaDetritus)(i, j, k, grid, ::Val{:N}, clock, fields)
    μᵖ = bgc.maximum_plankton_growth_rate
    μᵇ = bgc.maximum_bacteria_growth_rate
    r = bgc.linear_remineralization_rate
    gₘ = bgc.maximum_grazing_rate
    kᴰ = bgc.detritus_half_saturation
    kᴺ = bgc.nutrient_half_saturation
    kᵍ = bgc.grazing_half_saturation
    kᴵ = bgc.PAR_half_saturation
    λ = bgc.PAR_attenuation_scale
    y = bgc.bacteria_yield
    γ = bgc.zooplankton_yield

    # Available photosynthetic radiation
    z = znode(i, j, k, grid, c, c, c)

    # TODO: design a user interface for prescribing incoming shortwave
    I = 700 * exp(z / λ)

    P = @inbounds fields.P[i, j, k]
    Z = @inbounds fields.Z[i, j, k]
    D1 = @inbounds fields.D1[i, j, k] 
    D2 = @inbounds fields.D2[i, j, k]
    D = D1 .+ D2
    B = @inbounds fields.B[i, j, k]
    N = @inbounds fields.N[i, j, k]
    
    if sum(B) > 0
        return - phytoplankton_production(μᵖ, kᴺ, kᴵ, I, N, P) + bacteria_production(μᵇ, kᴰ, y, D, B) * (1/y - 1) 
               + zooplankton_graze_phytoplankton(gₘ, kᵍ, γ, P, Z) * (1/γ - 1) + zooplankton_graze_bacteria(gₘ, kᵍ, γ, B, Z) * (1/γ - 1)
    elseif sum(B) == 0
        return - phytoplankton_production(μᵖ, kᴺ, kᴵ, I, N, P) + detritus_remineralization(r, D)
        + zooplankton_graze_phytoplankton(gₘ, kᵍ, γ, P, Z) * (1/γ - 1)
    end

end

@inline function (bgc::NutrientsPlanktonBacteriaDetritus)(i, j, k, grid, ::Val{:P}, clock, fields)
    μᵖ = bgc.maximum_plankton_growth_rate
    gₘ = bgc.maximum_grazing_rate
    kᴺ = bgc.nutrient_half_saturation
    kᴵ = bgc.PAR_half_saturation
    kᵍ = bgc.grazing_half_saturation
    λ = bgc.PAR_attenuation_scale
    mlin = bgc.linear_mortality_rate
    mq = bgc.quadratic_mortality_rate
    γ = bgc.zooplankton_yield

    # Available photosynthetic radiation
    z = znode(i, j, k, grid, c, c, c)

    # TODO: design a user interface for prescribing incoming shortwave
    I = 700 * exp(z / λ)

    P = @inbounds fields.P[i, j, k]
    Z = @inbounds fields.Z[i, j, k]
    N = @inbounds fields.N[i, j, k]

    return phytoplankton_production(μᵖ, kᴺ, kᴵ, I, N, P) - phytoplankton_mortality(mlin, mq, P) - zooplankton_graze_phytoplankton(gₘ, kᵍ, γ, P, Z) / γ 
end

@inline function (bgc::NutrientsPlanktonBacteriaDetritus)(i, j, k, grid, ::Val{:Z}, clock, fields)
    gₘ = bgc.maximum_grazing_rate
    kᵍ = bgc.grazing_half_saturation
    mlin = bgc.linear_mortality_rate
    mq_Z = bgc.quadratic_mortality_rate_Z
    γ = bgc.zooplankton_yield

    P = @inbounds fields.P[i, j, k]
    B = @inbounds fields.B[i, j, k]
    Z = @inbounds fields.Z[i, j, k]

    return zooplankton_graze_phytoplankton(gₘ, kᵍ, γ, P, Z) + zooplankton_graze_bacteria(gₘ, kᵍ, γ, B, Z) - zooplankton_mortality(mlin, mq_Z, Z)   
end

@inline function (bgc::NutrientsPlanktonBacteriaDetritus)(i, j, k, grid, ::Val{:B}, clock, fields)
    μᵇ = bgc.maximum_bacteria_growth_rate
    gₘ = bgc.maximum_grazing_rate
    kᴰ = bgc.detritus_half_saturation
    kᵍ  = bgc.grazing_half_saturation
    mlin = bgc.linear_mortality_rate
    mq = bgc.quadratic_mortality_rate
    y = bgc.bacteria_yield
    γ = bgc.zooplankton_yield

    D1 = @inbounds fields.D1[i, j, k]
    D2 = @inbounds fields.D2[i, j, k]
    D = D1 .+ D2
    B = @inbounds fields.B[i, j, k]
    Z = @inbounds fields.Z[i, j, k]

    return bacteria_production(μᵇ, kᴰ, y, D, B) - bacteria_mortality(mlin, mq, B) - zooplankton_graze_bacteria(gₘ, kᵍ, γ, B, Z) / γ
end

@inline function (bgc::NutrientsPlanktonBacteriaDetritus)(i, j, k, grid, ::Val{:D1}, clock, fields)
    μᵇ = bgc.maximum_bacteria_growth_rate
    kᴰ = bgc.detritus_half_saturation
    r = bgc.linear_remineralization_rate
    y = bgc.bacteria_yield
    mlin = bgc.linear_mortality_rate
    mq = bgc.quadratic_mortality_rate
    mq_Z = bgc.quadratic_mortality_rate_Z

    P = @inbounds fields.P[i, j, k]
    Z = @inbounds fields.Z[i, j, k]
    D = @inbounds fields.D1[i, j, k] 
    B = @inbounds fields.B[i, j, k]

    if sum(B) > 0
        return bacteria_mortality(mlin, mq, B) + phytoplankton_mortality(mlin, mq, P) + zooplankton_mortality(mlin, mq_Z, Z) - bacteria_production(μᵇ, kᴰ, y, D, B) / y 
    elseif sum(B) == 0
        return phytoplankton_mortality(mlin, mq, P) + zooplankton_mortality(mlin, mq_Z, Z) - detritus_remineralization(r, D)
    end
end

@inline function (bgc::NutrientsPlanktonBacteriaDetritus)(i, j, k, grid, ::Val{:D2}, clock, fields)
    μᵇ = bgc.maximum_bacteria_growth_rate
    kᴰ = bgc.detritus_half_saturation
    r = bgc.linear_remineralization_rate
    y = bgc.bacteria_yield
    mlin = bgc.linear_mortality_rate
    mq = bgc.quadratic_mortality_rate
    mq_Z = bgc.quadratic_mortality_rate_Z

    P = @inbounds fields.P[i, j, k]
    Z = @inbounds fields.Z[i, j, k]
    D = @inbounds fields.D2[i, j, k]
    B = @inbounds fields.B[i, j, k]

    if sum(B) > 0
        return bacteria_mortality(mlin, mq, B) + phytoplankton_mortality(mlin, mq, P) + zooplankton_mortality(mlin, mq_Z, Z) - bacteria_production(μᵇ, kᴰ, y, D, B) / y
    elseif sum(B) == 0
        return phytoplankton_mortality(mlin, mq, P) + zooplankton_mortality(mlin, mq_Z, Z) - detritus_remineralization(r, D)
    end
end
