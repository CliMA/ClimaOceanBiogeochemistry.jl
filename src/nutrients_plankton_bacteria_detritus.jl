using Oceananigans.Units: day
using Oceananigans.Grids: znode, Center
using Oceananigans.Fields: ZeroField, ConstantField
using Oceananigans.Biogeochemistry: AbstractBiogeochemistry

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers, biogeochemical_drift_velocity

const c = Center()

"""
    NutrientsPlanktonBacteriaDetritus(; kw...)

Return a four-tracer biogeochemistry model for the interaction of nutrients, plankton
bacteria, and detritus.

Parameters
==========
  * `maximum_plankton_growth_rate`: (s⁻¹) Growth rate of plankton `P` unlimited by the
                                     availability of nutrients and light. Default: 1/day.

  * `maximum_bacteria_growth_rate`: (s⁻¹) Growth rate of plankton `B` unlimited by the
                                    availability of nutrients and light. Default = 0.5/day.

  * `bacteria_yield`: Determines fractional nutrient production by bacteria production
                      relative to consumption of detritus such that ``∂_t N / ∂_t D = 1 - y``,
                      where `y = bacteria_yield`. Default: 0.2.

  * `linear_mortality_rate`: (s⁻¹) Linear term of the mortality rate of both plankton and bacteria.

  * `quadratic_mortality_rate`: (s⁻¹) Quadratic term of the mortality rate of both plankton and bacteria.

  * `nutrient_half_saturation`: (mmol m⁻³) Half-saturation of nutrients for plankton production.

  * `detritus_half_saturation`: (mmol m⁻³) Half-saturation of nutrients for bacteria production.
                                Default = 10.0 mmol m⁻³.

  * `PAR_half_saturation`: (W m⁻²) Half-saturation of photosynthetically available radiation (PAR)
                           for plankton production.

  * `PAR_attenuation_scale`: (m) Depth scale over which photosynthetically available radiation (PAR)
                             attenuates exponentially.

  * `detritus_sinking_speed`: (m s⁻¹) Sinking velocity of detritus.

Tracer names
============
  * `N`: nutrients
  * `P`: plankton
  * `B`: bacteria
  * `D`: detritus

Biogeochemical functions
========================
  * transitions for `N`, `P`, `B`, `D`
  * `biogeochemical_drift_velocity` for `D`, modeling the sinking of detritus at
    a constant `detritus_sinking_speed`.
"""
Base.@kwdef struct NutrientsPlanktonBacteriaDetritus{FT} <: AbstractBiogeochemistry
    maximum_plankton_growth_rate :: FT = 1/day
    maximum_bacteria_growth_rate :: FT = 0.5/day
    bacteria_yield :: FT               = 0.2
    linear_mortality_rate :: FT        = 0.01/day # m³/mmol/day
    quadratic_mortality_rate :: FT     = 1/day # m³/mmol/day
    nutrient_half_saturation :: FT     = 0.1   # mmol m⁻³
    detritus_half_saturation :: FT     = 0.1   # mmol m⁻³
    PAR_half_saturation :: FT          = 10.0  # W m⁻²
    PAR_attenuation_scale :: FT        = 25.0  # m
    detritus_sinking_speed :: FT       = 10/day # m s⁻¹
end

const NPBD = NutrientsPlanktonBacteriaDetritus

@inline required_biogeochemical_tracers(::NPBD) = (:N, :P, :B, :D)

@inline function biogeochemical_drift_velocity(bgc::NPBD, ::Val{:D})
    u = ZeroField()
    v = ZeroField()
    w = ConstantField(- bgc.detritus_sinking_speed)
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

# Temperature modification to metabolic rates, following the Arrhenius equation
@inline temp_fun(Temp) = 0.8 .* exp.(-4000 .*(1 ./ (Temp .+ 273.15) .- 1 ./ 293.15))

@inline bacteria_production(μᵇ, kᴰ, D, B) = (μᵇ * D / (D + kᴰ) * B) * temp_fun(Temp)
@inline plankton_production(μᵖ, kᴺ, kᴵ, I, N, P) = (μᵖ * min(N / (N + kᴺ) , I / (I + kᴵ)) * P) * temp_fun(Temp)
@inline bacteria_mortality(mlin,mq, B) = mlin * B + mq * B^2
@inline plankton_mortality(mlin,mq, P) = mlin * P + mq * P^2

@inline function (bgc::NutrientsPlanktonBacteriaDetritus)(i, j, k, grid, ::Val{:N}, clock, fields)
    μᵖ = bgc.maximum_plankton_growth_rate
    μᵇ = bgc.maximum_bacteria_growth_rate
    kᴰ = bgc.detritus_half_saturation
    kᴺ = bgc.nutrient_half_saturation
    kᴵ = bgc.PAR_half_saturation
    λ = bgc.PAR_attenuation_scale
    mlin = bgc.linear_mortality_rate
    mq = bgc.quadratic_mortality_rate
    y = bgc.bacteria_yield

    # Available photosynthetic radiation
    z = znode(i, j, k, grid, c, c, c)

    # TODO: design a user interface for prescribing incoming shortwave
    I = 700 * exp(z / λ)

    P = @inbounds fields.P[i, j, k]
    D = @inbounds fields.D[i, j, k]
    B = @inbounds fields.B[i, j, k]
    N = @inbounds fields.N[i, j, k]
    
    return - plankton_production(μᵖ, kᴺ, kᴵ, I, N, P) + bacteria_production(μᵇ, kᴰ, D, B) * (1/y - 1)
end

@inline function (bgc::NutrientsPlanktonBacteriaDetritus)(i, j, k, grid, ::Val{:P}, clock, fields)
    μᵖ = bgc.maximum_plankton_growth_rate
    kᴺ = bgc.nutrient_half_saturation
    kᴵ = bgc.PAR_half_saturation
    λ = bgc.PAR_attenuation_scale
    mlin = bgc.linear_mortality_rate
    mq = bgc.quadratic_mortality_rate

    # Available photosynthetic radiation
    z = znode(i, j, k, grid, c, c, c)

    # TODO: design a user interface for prescribing incoming shortwave
    I = 700 * exp(z / λ)

    P = @inbounds fields.P[i, j, k]
    N = @inbounds fields.N[i, j, k]

    return plankton_production(μᵖ, kᴺ, kᴵ, I, N, P) - plankton_mortality(mlin, mq, P)
end

@inline function (bgc::NutrientsPlanktonBacteriaDetritus)(i, j, k, grid, ::Val{:B}, clock, fields)
    μᵇ = bgc.maximum_bacteria_growth_rate
    kᴰ = bgc.detritus_half_saturation
    mlin = bgc.linear_mortality_rate
    mq = bgc.quadratic_mortality_rate

    D = @inbounds fields.D[i, j, k]
    B = @inbounds fields.B[i, j, k]

    return bacteria_production(μᵇ, kᴰ, D, B) - bacteria_mortality(mlin, mq, B)
end

@inline function (bgc::NutrientsPlanktonBacteriaDetritus)(i, j, k, grid, ::Val{:D}, clock, fields)
    μᵇ = bgc.maximum_bacteria_growth_rate
    kᴰ = bgc.detritus_half_saturation
    y = bgc.bacteria_yield
    mlin = bgc.linear_mortality_rate
    mq = bgc.quadratic_mortality_rate

    P = @inbounds fields.P[i, j, k]
    D = @inbounds fields.D[i, j, k]
    B = @inbounds fields.B[i, j, k]

    return bacteria_mortality(mlin, mq, B) + plankton_mortality(mlin, mq, P) - bacteria_production(μᵇ, kᴰ, D, B) / y
end
