using Oceananigans
using Oceananigans.Architectures: arch_array, architecture
using Oceananigans.Units: day
using Oceananigans.Grids: znode, Center, AbstractTopology, Flat, Bounded
using Oceananigans.BoundaryConditions: ImpenetrableBoundaryCondition, fill_halo_regions!
using Oceananigans.Fields: ZeroField, ZFaceField
using Oceananigans.Biogeochemistry: AbstractBiogeochemistry

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers, biogeochemical_drift_velocity

const c = Center()

struct MultiNPZBD{P, B, Z, D, BDC, PNC, ZPC, FT, W} <: AbstractBiogeochemistry
    nutrient_types :: Int
    plankton_populations :: Int
    zooplankton_populations :: Int
    bacteria_populations :: Int
    detritus_types :: Int
    maximum_plankton_growth_rate :: P
    maximum_bacteria_growth_rate :: B
    maximum_grazing_rate :: Z
    bacteria_yield :: B
    zooplankton_yield :: Z
    linear_remineralization_rate :: D
    linear_plankton_mortality_rate :: P
    linear_bacteria_mortality_rate :: B
    linear_zooplankton_mortality_rate :: Z
    quadratic_plankton_mortality_rate :: P
    quadratic_bacteria_mortality_rate :: B
    quadratic_zooplankton_mortality_rate :: Z
    nutrient_half_saturation :: P
    detritus_half_saturation :: B
    grazing_half_saturation  :: Z
    bacteria_detritus_connectivity :: BDC
    plankton_nutrient_connectivity :: PNC
    zooplankton_plankton_connectivity :: ZPC
    PAR_half_saturation :: FT
    PAR_attenuation_scale :: FT
    detritus_vertical_velocity :: W
end

arrayify(grid, n, t::Number) = arch_array(architecture(grid), [convert(eltype(grid),t) for _ = 1:n])
arrayify(grid, n, t::Array)  = arch_array(architecture(grid), t)
arrayify(grid, n, t)         = arch_array(architecture(grid), [t for _ = 1:n])

"""
    MultiNPZBD(grid; kw...)

Return a five-tracer biogeochemistry model for the interaction of nutrients (N), phytoplankton (P), 
zooplankton(Z), bacteria (B), and detritus (D).

Parameters
==========
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

* `detritus_sinking_speed`: (m s⁻¹) Sinking velocity of detritus.

Tracer names
============
* `N`: nutrients
* `P`: phytoplankton
* `Z`: zooplankton
* `B`: bacteria
* `D`: detritus 

Biogeochemical functions
========================
* transitions for `N`, `P`, `Z`, `B`, `D`
* `biogeochemical_drift_velocity` for `D`, modeling the sinking of detritus at
  a constant `detritus_sinking_speed`.

Example
=======

```julia
bgc = MultiNPZBD(grid, plankton_populations=2)
bgc.maximum_plankton_growth_rate[1] = 1/day
bgc.maximum_plankton_growth_rate[2] = 2/day
```
"""
function MultiNPZBD(grid;
                    nutrient_types = 1,
                    plankton_populations = 1,
                    bacteria_populations = 1,
                    zooplankton_populations = 1,
                    detritus_types  = 1,
                    maximum_plankton_growth_rate         = 1/day, # Add reference for each parameter
                    maximum_bacteria_growth_rate         = 1/day,
                    maximum_grazing_rate                 = 3/day,
                    bacteria_yield                       = 0.2,
                    zooplankton_yield                    = 0.3,
                    linear_remineralization_rate         = 0.03/day, 
                    linear_plankton_mortality_rate       = 0.01/day, # m³/mmol/day
                    linear_bacteria_mortality_rate       = 0.01/day,
                    linear_zooplankton_mortality_rate    = 0.01/day, 
                    quadratic_plankton_mortality_rate    = 0.1/day, # m³/mmol/day (zooplankton quadratic mortality)= 0.1,   # mmol m⁻³
                    quadratic_bacteria_mortality_rate    = 0.1/day, # mmol m⁻³
                    quadratic_zooplankton_mortality_rate = 1/day,   # mmol m⁻³
                    nutrient_half_saturation             = 0.1,
                    detritus_half_saturation             = 0.1,
                    grazing_half_saturation              = 3.0,
                    bacteria_detritus_connectivity       = nothing,
                    plankton_nutrient_connectivity       = nothing,
                    zooplankton_plankton_connectivity    = nothing,
                    PAR_half_saturation                  = 10.0,
                    PAR_attenuation_scale                = 25.0,
                    detritus_vertical_velocity           = -10/day) # m s⁻¹

    FT = eltype(grid)

    NN = nutrient_types
    NP = plankton_populations
    NB = bacteria_populations
    NZ = zooplankton_populations
    ND = detritus_types

    maximum_plankton_growth_rate         = arrayify(grid, NP, maximum_plankton_growth_rate)
    maximum_bacteria_growth_rate         = arrayify(grid, NB, maximum_bacteria_growth_rate)
    maximum_grazing_rate                 = arrayify(grid, NZ, maximum_grazing_rate)
    bacteria_yield                       = arrayify(grid, NB, bacteria_yield)
    zooplankton_yield                    = arrayify(grid, NZ, zooplankton_yield)
    linear_remineralization_rate         = arrayify(grid, ND, linear_remineralization_rate)
    linear_plankton_mortality_rate       = arrayify(grid, NP, linear_plankton_mortality_rate)
    linear_bacteria_mortality_rate       = arrayify(grid, NB, linear_bacteria_mortality_rate)
    linear_zooplankton_mortality_rate    = arrayify(grid, NZ, linear_zooplankton_mortality_rate)
    quadratic_plankton_mortality_rate    = arrayify(grid, NP, quadratic_plankton_mortality_rate)
    quadratic_bacteria_mortality_rate    = arrayify(grid, NB, quadratic_bacteria_mortality_rate)
    quadratic_zooplankton_mortality_rate = arrayify(grid, NZ, quadratic_zooplankton_mortality_rate)
    nutrient_half_saturation             = arrayify(grid, NP, nutrient_half_saturation)
    detritus_half_saturation             = arrayify(grid, NB, detritus_half_saturation)
    grazing_half_saturation              = arrayify(grid, NZ, grazing_half_saturation)

    if detritus_vertical_velocity isa Number
        w₀ = detritus_vertical_velocity 
        no_penetration = ImpenetrableBoundaryCondition()

        bcs = FieldBoundaryConditions(grid, (Center, Center, Face),
                                      top=no_penetration, bottom=no_penetration)

        detritus_vertical_velocity = ZFaceField(grid, boundary_conditions = bcs)

        set!(detritus_vertical_velocity, w₀)

        fill_halo_regions!(detritus_vertical_velocity)
    end

    detritus_vertical_velocity = arrayify(grid, NZ, detritus_vertical_velocity)

    FT = eltype(grid)

    npzbd = MultiNPZBD(NN, NP, NZ, NB, ND,
                       maximum_plankton_growth_rate,
                       maximum_bacteria_growth_rate,
                       maximum_grazing_rate,
                       bacteria_yield,
                       zooplankton_yield,
                       linear_remineralization_rate,
                       linear_plankton_mortality_rate,
                       linear_bacteria_mortality_rate,
                       linear_zooplankton_mortality_rate,
                       quadratic_plankton_mortality_rate,
                       quadratic_bacteria_mortality_rate,
                       quadratic_zooplankton_mortality_rate,
                       nutrient_half_saturation,
                       detritus_half_saturation,
                       grazing_half_saturation,
                       bacteria_detritus_connectivity,
                       plankton_nutrient_connectivity,
                       zooplankton_plankton_connectivity,
                       convert(FT, PAR_half_saturation),
                       convert(FT, PAR_attenuation_scale),
                       detritus_vertical_velocity)

    return npzbd
end

# Add tracer names if there are multiple groups in each tracer (e.g. D1, D2)
@inline function required_biogeochemical_tracers(bgc::MultiNPZBD)
    NN = bgc.nutrient_types
    NP = bgc.plankton_populations
    NZ = bgc.zooplankton_populations
    NB = bgc.bacteria_populations
    ND = bgc.detritus_types

    nutrient_names    = Tuple(Symbol(:N, n) for n = 1:NN)
    plankton_names    = Tuple(Symbol(:P, n) for n = 1:NP)
    zooplankton_names = Tuple(Symbol(:Z, n) for n = 1:NZ)
    bacteria_names    = Tuple(Symbol(:B, n) for n = 1:NB)
    detritus_names    = Tuple(Symbol(:D, n) for n = 1:ND)

    return tuple(nutrient_names...,
                 plankton_names...,
                 zooplankton_names...,
                 bacteria_names...,
                 detritus_names...)   
end

# @inline required_biogeochemical_tracers(::MultiNPZBD) = required_biogeochemical_tracers(bgc)

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

@inline function phytoplankton_production(μᵖ, kᴺ, kᴵ, I, N, P)
    Ñ = N / (N + kᴺ)
    Ĩ = I / (I + kᴵ)
    return μᵖ * P * min(Ñ, Ĩ)
end

@inline zooplankton_graze_phytoplankton(gₘ, kᵍ, γ, P, Z) = γ * gₘ * P / (P + kᵍ) * Z
@inline zooplankton_graze_bacteria(gₘ, kᵍ, γ, B, Z) = γ * gₘ * B / (B + kᵍ) * Z
@inline bacteria_mortality(mlinB, mqB, B) = mlinB * B + mqB * B^2
@inline phytoplankton_mortality(mlinP, mqP, P) = mlinP * P + mqP * P^2
@inline zooplankton_mortality(mlinZ, mqZ, Z) = mlinZ * Z + mqZ * Z^2
@inline detritus_remineralization(r, D) = r * D


for C in (:N, :P, :Z, :B, :D)
    @eval begin
        @inline numbered_tracer_tuple(::Val{$(Meta.quot(C))}, ::Val{1}, fields) = tuple(getproperty(fields, Symbol($(Meta.quot(C)), 1)))

        @inline numbered_tracer_tuple(::Val{$(Meta.quot(C))}, ::Val{2}, fields) = (getproperty(fields, Symbol($(Meta.quot(C)), 1)),
                                                                                   getproperty(fields, Symbol($(Meta.quot(C)), 2)))

        @inline numbered_tracer_tuple(::Val{$(Meta.quot(C))}, ::Val{3}, fields) = (getproperty(fields, Symbol($(Meta.quot(C)), 1)),
                                                                                   getproperty(fields, Symbol($(Meta.quot(C)), 2)),
                                                                                   getproperty(fields, Symbol($(Meta.quot(C)), 3)))

        @inline numbered_tracer_tuple(::Val{$(Meta.quot(C))}, ::Val{4}, fields) = (getproperty(fields, Symbol($(Meta.quot(C)), 1)),
                                                                                   getproperty(fields, Symbol($(Meta.quot(C)), 2)),
                                                                                   getproperty(fields, Symbol($(Meta.quot(C)), 3)),
                                                                                   getproperty(fields, Symbol($(Meta.quot(C)), 4)))

        @inline numbered_tracer_tuple(::Val{$(Meta.quot(C))}, ::Val{5}, fields) = (getproperty(fields, Symbol($(Meta.quot(C)), 1)),
                                                                                   getproperty(fields, Symbol($(Meta.quot(C)), 2)),
                                                                                   getproperty(fields, Symbol($(Meta.quot(C)), 3)),
                                                                                   getproperty(fields, Symbol($(Meta.quot(C)), 4)),
                                                                                   getproperty(fields, Symbol($(Meta.quot(C)), 5)))

    end
end

@inline function tupleindex(i, j, k, t::NTuple{N}) where N
    ntuple(Val(N)) do n
        @inbounds t[n][i, j, k]
    end
end

# TODO: How to update ::Val{:Dn} as the input?
# Change this number and recompile if you would like more detritus types
MAX_TYPES = 5

for n = 1:MAX_TYPES
    @eval begin

        @inline function biogeochemical_drift_velocity(bgc::MultiNPZBD, ::Val{:D$n})
            u = ZeroField()
            v = ZeroField()
            w = bgc.detritus_vertical_velocity[$n]
            return (; u, v, w)
        end

        @inline function (bgc::MultiNPZBD)(i, j, k, grid, ::Val{:N$n}, clock, fields)
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

            Nᴺ = length(kᴺ)
            Nᴾ = length(μᵖ)
            Nᶻ = length(gₘ)
            Nᴰ = length(kᴰ)
            Nᴮ = length(μᵇ)

            N = numbered_tracer_tuple(Val(:N), Val(Nᴺ), fields)
            P = numbered_tracer_tuple(Val(:P), Val(Nᴾ), fields)
            Z = numbered_tracer_tuple(Val(:Z), Val(Nᶻ), fields)
            D = numbered_tracer_tuple(Val(:D), Val(Nᴰ), fields)
            B = numbered_tracer_tuple(Val(:B), Val(Nᴮ), fields)

            Nⁱʲᵏ = tupleindex(N)
            Pⁱʲᵏ = tupleindex(P)
            Zⁱʲᵏ = tupleindex(Z)
            Dⁱʲᵏ = tupleindex(D)
            Bⁱʲᵏ = tupleindex(B)

            # μᵖ

            if sum(B) > 0
                return (- phytoplankton_production(μᵖ, kᴺ, kᴵ, I, N, P)
                        .+ bacteria_production(μᵇ, kᴰ, y, D, B) * (1 / y - 1)
                        .+ zooplankton_graze_phytoplankton(gₘ, kᵍ, γ, P, Z) * (1 / γ - 1)
                        .+ zooplankton_graze_bacteria(gₘ, kᵍ, γ, B, Z) * (1 / γ - 1))

            elseif sum(B) == 0
                return (- phytoplankton_production(μᵖ, kᴺ, kᴵ, I, N, P)
                        .+ detritus_remineralization(r, D)
                        .+ zooplankton_graze_phytoplankton(gₘ, kᵍ, γ, P, Z) * (1 / γ - 1))
            end 
        end

        @inline function (bgc::MultiNPZBD)(i, j, k, grid, ::Val{:P$n}, clock, fields)
            μᵖ = bgc.maximum_plankton_growth_rate
            gₘ = bgc.maximum_grazing_rate
            kᴺ = bgc.nutrient_half_saturation
            kᴵ = bgc.PAR_half_saturation
            kᵍ = bgc.grazing_half_saturation
            λ = bgc.PAR_attenuation_scale
            mlinP = bgc.linear_plankton_mortality_rate
            mqP = bgc.quadratic_plankton_mortality_rate
            γ = bgc.zooplankton_yield

            # Available photosynthetic radiation
            z = znode(i, j, k, grid, c, c, c)

            # TODO: design a user interface for prescribing incoming shortwave
            I = 700 * exp(z / λ)

            Nᴾ = length(μᵖ)

            P = @inbounds fields.P1[i, j, k]
            Z = @inbounds fields.Z1[i, j, k]
            N = @inbounds fields.N1[i, j, k]

            return (phytoplankton_production(μᵖ, kᴺ, kᴵ, I, N, P)
                    - phytoplankton_mortality(mlinP, mqP, P)
                    - zooplankton_graze_phytoplankton(gₘ, kᵍ, γ, P, Z) / γ)
        end

        @inline function (bgc::MultiNPZBD)(i, j, k, grid, ::Val{:Z$n}, clock, fields)
            gₘ = bgc.maximum_grazing_rate
            kᵍ = bgc.grazing_half_saturation
            mlinZ = bgc.linear_zooplankton_mortality_rate
            mqZ = bgc.quadratic_zooplankton_mortality_rate
            γ = bgc.zooplankton_yield

            P = @inbounds fields.P1[i, j, k]
            B = @inbounds fields.B1[i, j, k]
            Z = @inbounds fields.Z1[i, j, k]

            return (zooplankton_graze_phytoplankton(gₘ, kᵍ, γ, P, Z)
                    + zooplankton_graze_bacteria(gₘ, kᵍ, γ, B, Z)
                    - zooplankton_mortality(mlinZ, mqZ, Z))
        end

        @inline function (bgc::MultiNPZBD)(i, j, k, grid, ::Val{:B$n}, clock, fields)
            μᵇ = bgc.maximum_bacteria_growth_rate
            gₘ = bgc.maximum_grazing_rate
            kᴰ = bgc.detritus_half_saturation
            kᵍ  = bgc.grazing_half_saturation
            mlinB = bgc.linear_bacteria_mortality_rate
            mqB = bgc.quadratic_bacteria_mortality_rate
            y = bgc.bacteria_yield
            γ = bgc.zooplankton_yield

            D = @inbounds fields.D1[i, j, k] 
            B = @inbounds fields.B1[i, j, k]
            Z = @inbounds fields.Z1[i, j, k]

            return (bacteria_production(μᵇ, kᴰ, y, D, B)
                    - bacteria_mortality(mlinB, mqB, B)
                    - zooplankton_graze_bacteria(gₘ, kᵍ, γ, B, Z) / γ) 
        end

        @inline function (bgc::MultiNPZBD)(i, j, k, grid, ::Val{:D$n}, clock, fields)
            μᵇ = bgc.maximum_bacteria_growth_rate
            kᴰ = bgc.detritus_half_saturation
            r = bgc.linear_remineralization_rate
            y = bgc.bacteria_yield
            mlinP = bgc.linear_plankton_mortality_rate
            mqP = bgc.quadratic_plankton_mortality_rate
            mlinZ = bgc.linear_zooplankton_mortality_rate
            mqZ = bgc.quadratic_zooplankton_mortality_rate
            mlinB = bgc.linear_bacteria_mortality_rate
            mqB = bgc.quadratic_bacteria_mortality_rate

            P = @inbounds fields.P1[i, j, k]
            Z = @inbounds fields.Z1[i, j, k]
            D = @inbounds fields.D1[i, j, k] 
            B = @inbounds fields.B1[i, j, k]

            if sum(B) > 0
                return (bacteria_mortality(mlinB, mqB, B)
                        + phytoplankton_mortality(mlinP, mqP, P)
                        + zooplankton_mortality(mlinZ, mqZ, Z)
                        - bacteria_production(μᵇ, kᴰ, y, D, B) / y)
            elseif sum(B) == 0
                return (phytoplankton_mortality(mlinP, mqP, P)
                        + zooplankton_mortality(mlinZ, mqZ, Z)
                        - detritus_remineralization(r, D))
            end 
        end
    end
end

#=
@inline function (bgc::MultiNPZBD)(i, j, k, grid, ::Val{:D2}, clock, fields)
    μᵇ = bgc.maximum_bacteria_growth_rate
    kᴰ = bgc.detritus_half_saturation
    r = bgc.linear_remineralization_rate
    y = bgc.bacteria_yield
    mlinP = bgc.linear_plankton_mortality_rate
    mqP = bgc.quadratic_plankton_mortality_rate
    mlinZ = bgc.linear_zooplankton_mortality_rate
    mqZ = bgc.quadratic_zooplankton_mortality_rate
    mlinB = bgc.linear_bacteria_mortality_rate
    mqB = bgc.quadratic_bacteria_mortality_rate

    P = @inbounds fields.P1[i, j, k]
    Z = @inbounds fields.Z1[i, j, k]
    D = @inbounds fields.D2[i, j, k] 
    B = @inbounds fields.B1[i, j, k]

    if sum(B) > 0
        return bacteria_mortality(mlinB, mqB, B) + phytoplankton_mortality(mlinP, mqP, P) + zooplankton_mortality(mlinZ, mqZ, Z) - bacteria_production(μᵇ, kᴰ, y, D, B) / y 
    elseif sum(B) == 0
        return phytoplankton_mortality(mlinP, mqP, P) + zooplankton_mortality(mlinZ, mqZ, Z) - detritus_remineralization(r, D)
    end 
end
=#
