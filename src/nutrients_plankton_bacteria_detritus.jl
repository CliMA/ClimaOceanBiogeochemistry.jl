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

Tracer names
============
  - `N`: nutrients
  - `P`: plankton
  - `B`: bacteria
  - `D`: detritus

Biogeochemical functions
========================
  - transitions for `N`, `P`, `B`, `D`
  - `biogeochemical_drift_velocity` for `D`, modeling the sinking of detritus at
    a constant `detritus_sinking_speed`.
"""
Base.@kwdef struct NutrientsPlanktonBacteriaDetritus{FT} <: AbstractBiogeochemistry
    maximum_plankton_growth_rate :: FT = 1/day
    maximum_bacteria_growth_rate :: FT = 0.5/day
    bacteria_yield :: FT               = 0.2
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

@inline bacteria_production(μᵇ, kᴰ, D, B) = μᵇ * D / (D + kᴰ) * B
@inline plankton_production(μᵖ, kᴺ, kᴵ, I, N, P) = μᵖ * (N / (N + kᴺ) * I / (I + kᴵ)) * P
@inline bacteria_mortality(m, B) = m * B^2
@inline plankton_mortality(m, P) = m * P^2

@inline function (bgc::NutrientsPlanktonBacteriaDetritus)(i, j, k, grid, ::Val{:N}, clock, fields)
    μᵖ = bgc.maximum_plankton_growth_rate
    μᵇ = bgc.maximum_bacteria_growth_rate
    kᴰ = bgc.detritus_half_saturation
    kᴺ = bgc.nutrient_half_saturation
    kᴵ = bgc.PAR_half_saturation
    λ = bgc.PAR_attenuation_scale
    m = bgc.quadratic_mortality_rate
    y = bgc.bacteria_yield

    # Available photosynthetic radiation
    z = znode(i, j, k, grid, c, c, c)
    I = exp(z / λ)

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
    m = bgc.quadratic_mortality_rate

    # Available photosynthetic radiation
    z = znode(i, j, k, grid, c, c, c)
    I = exp(z / λ)

    P = @inbounds fields.P[i, j, k]
    N = @inbounds fields.N[i, j, k]

    return plankton_production(μᵖ, kᴺ, kᴵ, I, N, P) - plankton_mortality(m, P)
end

@inline function (bgc::NutrientsPlanktonBacteriaDetritus)(i, j, k, grid, ::Val{:B}, clock, fields)
    μᵇ = bgc.maximum_bacteria_growth_rate
    kᴰ = bgc.detritus_half_saturation
    m = bgc.quadratic_mortality_rate

    D = @inbounds fields.D[i, j, k]
    B = @inbounds fields.B[i, j, k]

    return bacteria_production(μᵇ, kᴰ, D, B) - bacteria_mortality(m, B)
end

@inline function (bgc::NutrientsPlanktonBacteriaDetritus)(i, j, k, grid, ::Val{:D}, clock, fields)
    μᵇ = bgc.maximum_bacteria_growth_rate
    kᴰ = bgc.detritus_half_saturation
    y = bgc.bacteria_yield
    m = bgc.quadratic_mortality_rate

    P = @inbounds fields.P[i, j, k]
    D = @inbounds fields.D[i, j, k]
    B = @inbounds fields.B[i, j, k]

    return bacteria_mortality(m, B) + plankton_mortality(m, P) - bacteria_production(μᵇ, kᴰ, D, B) / y
end

