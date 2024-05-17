using Oceananigans
#using Oceananigans.Architectures: arch_array, architecture, on_architecture
using Oceananigans.Units: day
using Oceananigans.Grids: znode, Center, AbstractTopology, Flat, Bounded
using Oceananigans.BoundaryConditions: ImpenetrableBoundaryCondition, fill_halo_regions!
using Oceananigans.Fields: ZeroField, ZFaceField
using Oceananigans.Biogeochemistry: AbstractBiogeochemistry

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers, biogeochemical_drift_velocity

const c = Center()

"""
    MultiNPZBD(; kw...)

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
"""
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

"""
    MultiNPZBD(grid)

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
                       PAR_half_saturation,            
                       PAR_attenuation_scale,          
                       detritus_vertical_velocity)

    return npzbd
end

# Add tracer names if there are multiple groups in each tracer (e.g. D1, D2)
#=
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
=#
@inline required_biogeochemical_tracers(::MultiNPZBD) = (:N1,:P1,:P2,:Z1,:Z2,:B1,:B2,:D1,:D2,:D3,:D4,:D5)

# Set sinking velocities of POM
@inline function biogeochemical_drift_velocity(bgc::MultiNPZBD, ::Val{:D4})
    u = ZeroField()
    v = ZeroField()
    w = bgc.detritus_vertical_velocity./2
    return (; u, v, w)
end
@inline function biogeochemical_drift_velocity(bgc::MultiNPZBD, ::Val{:D5})
    u = ZeroField()
    v = ZeroField()
    w = bgc.detritus_vertical_velocity.*5
    return (; u, v, w)
end

@inline bacteria_production(μᵇ, kᴰ, y, D, B) = y .* μᵇ .* D ./ (D .+ kᴰ) .* B 
@inline phytoplankton_production(μᵖ, kᴺ, kᴵ, I, N, P) = (μᵖ * min(N / (N .+ kᴺ) , I / (I .+ kᴵ)) * P) 
@inline zooplankton_graze_phytoplankton(gₘ, kᵍ, γ, P, Z) = γ .* gₘ .* P ./ (P .+ kᵍ) .* Z
@inline zooplankton_graze_bacteria(gₘ, kᵍ, γ, B, Z) = γ .* gₘ .* B ./ (B .+ kᵍ) .* Z
@inline bacteria_mortality(mlinB,mqB, B) = mlinB .* B .+ mqB .* B^2
@inline phytoplankton_mortality(mlinP,mqP, P) = mlinP .* P .+ mqP .* P^2
@inline zooplankton_mortality(mlinZ,mqZ, Z) = mlinZ .* Z .+ mqZ .* Z^2
@inline detritus_remineralization(r, D) = r .* D

@inline function (bgc::MultiNPZBD)(i, j, k, grid, ::Val{:N1}, clock, fields)
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

    N1 = @inbounds fields.N1[i, j, k]
    P1 = @inbounds fields.P1[i, j, k] 
    P2 = @inbounds fields.P2[i, j, k] 
    Z1 = @inbounds fields.Z1[i, j, k]
    Z2 = @inbounds fields.Z2[i, j, k]
    B1 = @inbounds fields.B1[i, j, k]
    B2 = @inbounds fields.B2[i, j, k]
    D1 = @inbounds fields.D1[i, j, k] 
    D2 = @inbounds fields.D2[i, j, k] 
    D3 = @inbounds fields.D3[i, j, k] 
    D4 = @inbounds fields.D4[i, j, k] 
    D5 = @inbounds fields.D5[i, j, k] 
    
    if sum(B1+B2) > 0
        return (.- phytoplankton_production(μᵖ, kᴺ, kᴵ, I, N1, P1).- phytoplankton_production(μᵖ*2, kᴺ, kᴵ, I, N1, P2) .+ bacteria_production(μᵇ*2, kᴰ, y, D1, B1) * (1/y .- 1) .+ bacteria_production(μᵇ, kᴰ, y, D4, B1) * (1/y .- 1) .+ bacteria_production(μᵇ*2, kᴰ, y, D2, B2) * (1/y .- 1) .+ bacteria_production(μᵇ*2/2, kᴰ, y, D3, B2) * (1/y .- 1) .+ bacteria_production(μᵇ*2, kᴰ, y, D5, B2) * (1/y .- 1) .+ zooplankton_graze_phytoplankton(gₘ, kᵍ, γ, P1, Z1) * (1/γ .- 1) .+ zooplankton_graze_phytoplankton(gₘ, kᵍ, γ, P2, Z1) * (1/γ .- 1) .+ zooplankton_graze_bacteria(gₘ, kᵍ, γ, B1, Z2) * (1/γ .- 1) .+ zooplankton_graze_bacteria(gₘ, kᵍ, γ, B2, Z2) * (1/γ .- 1) )
    elseif sum(B1+B2) == 0
        return (.- phytoplankton_production(μᵖ, kᴺ, kᴵ, I, N1, P1).- phytoplankton_production(μᵖ*2, kᴺ, kᴵ, I, N1, P2) .+ detritus_remineralization(r*2, D1) .+detritus_remineralization(r, D2) .+detritus_remineralization(r/2, D3) .+detritus_remineralization(r, D4) .+detritus_remineralization(r, D5) .+zooplankton_graze_phytoplankton(gₘ, kᵍ, γ, P1, Z1) * (1/γ .- 1) .+ zooplankton_graze_phytoplankton(gₘ, kᵍ, γ, P2, Z1) * (1/γ .- 1))
    end
end

@inline function (bgc::MultiNPZBD)(i, j, k, grid, ::Val{:P1}, clock, fields)
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

    P1 = @inbounds fields.P1[i, j, k]
    Z1 = @inbounds fields.Z1[i, j, k]
    N1 = @inbounds fields.N1[i, j, k]

    return phytoplankton_production(μᵖ, kᴺ, kᴵ, I, N1, P1) - phytoplankton_mortality(mlinP, mqP, P1) - zooplankton_graze_phytoplankton(gₘ, kᵍ, γ, P1, Z1) /γ 
end

@inline function (bgc::MultiNPZBD)(i, j, k, grid, ::Val{:P2}, clock, fields)
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

    P2 = @inbounds fields.P2[i, j, k]
    Z1 = @inbounds fields.Z1[i, j, k]
    N1 = @inbounds fields.N1[i, j, k]

    return phytoplankton_production(μᵖ*2, kᴺ, kᴵ, I, N1, P2) - phytoplankton_mortality(mlinP, mqP, P2) - zooplankton_graze_phytoplankton(gₘ, kᵍ, γ, P2, Z1) /γ 
end

@inline function (bgc::MultiNPZBD)(i, j, k, grid, ::Val{:Z1}, clock, fields)
    gₘ = bgc.maximum_grazing_rate
    kᵍ = bgc.grazing_half_saturation
    mlinZ = bgc.linear_zooplankton_mortality_rate
    mqZ = bgc.quadratic_zooplankton_mortality_rate
    γ = bgc.zooplankton_yield

    P1 = @inbounds fields.P1[i, j, k] 
    P2 = @inbounds fields.P2[i, j, k]
    Z1 = @inbounds fields.Z1[i, j, k]

    return zooplankton_graze_phytoplankton(gₘ, kᵍ, γ, P1, Z1)+zooplankton_graze_phytoplankton(gₘ, kᵍ, γ, P2, Z1) - zooplankton_mortality(mlinZ, mqZ, Z1)   
end

@inline function (bgc::MultiNPZBD)(i, j, k, grid, ::Val{:Z2}, clock, fields)
    gₘ = bgc.maximum_grazing_rate
    kᵍ = bgc.grazing_half_saturation
    mlinZ = bgc.linear_zooplankton_mortality_rate
    mqZ = bgc.quadratic_zooplankton_mortality_rate
    γ = bgc.zooplankton_yield

    B1 = @inbounds fields.B1[i, j, k]
    B2 = @inbounds fields.B2[i, j, k]
    Z2 = @inbounds fields.Z2[i, j, k]

    return zooplankton_graze_bacteria(gₘ, kᵍ, γ, B1, Z2)+zooplankton_graze_bacteria(gₘ, kᵍ, γ, B2, Z2) - zooplankton_mortality(mlinZ, mqZ, Z2)   
end

@inline function (bgc::MultiNPZBD)(i, j, k, grid, ::Val{:B1}, clock, fields)
    μᵇ = bgc.maximum_bacteria_growth_rate
    gₘ = bgc.maximum_grazing_rate
    kᴰ = bgc.detritus_half_saturation
    kᵍ  = bgc.grazing_half_saturation
    mlinB = bgc.linear_bacteria_mortality_rate
    mqB = bgc.quadratic_bacteria_mortality_rate
    y = bgc.bacteria_yield
    γ = bgc.zooplankton_yield

    D1 = @inbounds fields.D1[i, j, k] 
    D4 = @inbounds fields.D4[i, j, k] 
    B1 = @inbounds fields.B1[i, j, k]
    Z2 = @inbounds fields.Z2[i, j, k]

    return (bacteria_production(μᵇ*2, kᴰ, y, D1, B1) + bacteria_production(μᵇ, kᴰ, y, D4, B1) - bacteria_mortality(mlinB, mqB, B1) - zooplankton_graze_bacteria(gₘ, kᵍ, γ, B1, Z2) /γ )
end

@inline function (bgc::MultiNPZBD)(i, j, k, grid, ::Val{:B2}, clock, fields)
    μᵇ = bgc.maximum_bacteria_growth_rate
    gₘ = bgc.maximum_grazing_rate
    kᴰ = bgc.detritus_half_saturation
    kᵍ  = bgc.grazing_half_saturation
    mlinB = bgc.linear_bacteria_mortality_rate
    mqB = bgc.quadratic_bacteria_mortality_rate
    y = bgc.bacteria_yield
    γ = bgc.zooplankton_yield

    D2 = @inbounds fields.D2[i, j, k] 
    D3 = @inbounds fields.D3[i, j, k] 
    D5 = @inbounds fields.D5[i, j, k] 
    B2 = @inbounds fields.B2[i, j, k]
    Z2 = @inbounds fields.Z2[i, j, k]

    return (bacteria_production(μᵇ*2, kᴰ, y, D2, B2) + bacteria_production(μᵇ*2/2, kᴰ, y, D3, B2) + bacteria_production(μᵇ*2, kᴰ, y, D5, B2) - bacteria_mortality(mlinB, mqB, B2) - zooplankton_graze_bacteria(gₘ, kᵍ, γ, B2, Z2) /γ )
end

@inline function (bgc::MultiNPZBD)(i, j, k, grid, ::Val{:D1}, clock, fields)
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

    P1 = @inbounds fields.P1[i, j, k] 
    P2 = @inbounds fields.P2[i, j, k] 
    Z1 = @inbounds fields.Z1[i, j, k] 
    Z2 = @inbounds fields.Z2[i, j, k]
    D1 = @inbounds fields.D1[i, j, k] 
    B1 = @inbounds fields.B1[i, j, k]
    B2 = @inbounds fields.B2[i, j, k]

    if sum(B1+B2) > 0
        return 0.2*(bacteria_mortality(mlinB, mqB, B1) + bacteria_mortality(mlinB, mqB, B2) + phytoplankton_mortality(mlinP, mqP, P1) + phytoplankton_mortality(mlinP, mqP, P2) + zooplankton_mortality(mlinZ, mqZ, Z1) + zooplankton_mortality(mlinZ, mqZ, Z2)) - bacteria_production(μᵇ*2, kᴰ, y, D1, B1) / y 
    elseif sum(B1+B2) == 0
        return 0.2*(phytoplankton_mortality(mlinP, mqP, P1) + phytoplankton_mortality(mlinP, mqP, P2) + zooplankton_mortality(mlinZ, mqZ, Z1) + zooplankton_mortality(mlinZ, mqZ, Z2)) - detritus_remineralization(r*2, D1)
    end 
end

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

    P1 = @inbounds fields.P1[i, j, k] 
    P2 = @inbounds fields.P2[i, j, k] 
    Z1 = @inbounds fields.Z1[i, j, k] 
    Z2 = @inbounds fields.Z2[i, j, k]
    D2 = @inbounds fields.D2[i, j, k] 
    B1 = @inbounds fields.B1[i, j, k]
    B2 = @inbounds fields.B2[i, j, k]

    if sum(B1+B2) > 0
        return 0.2*(bacteria_mortality(mlinB, mqB, B1) + bacteria_mortality(mlinB, mqB, B2) + phytoplankton_mortality(mlinP, mqP, P1) + phytoplankton_mortality(mlinP, mqP, P2) + zooplankton_mortality(mlinZ, mqZ, Z1) + zooplankton_mortality(mlinZ, mqZ, Z2)) - bacteria_production(μᵇ*2, kᴰ, y, D2, B2) / y 
    elseif sum(B1+B2) == 0
        return 0.2*(phytoplankton_mortality(mlinP, mqP, P1) + phytoplankton_mortality(mlinP, mqP, P2) + zooplankton_mortality(mlinZ, mqZ, Z1) + zooplankton_mortality(mlinZ, mqZ, Z2)) - detritus_remineralization(r, D2)
    end 
end

@inline function (bgc::MultiNPZBD)(i, j, k, grid, ::Val{:D3}, clock, fields)
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

    P1 = @inbounds fields.P1[i, j, k] 
    P2 = @inbounds fields.P2[i, j, k] 
    Z1 = @inbounds fields.Z1[i, j, k] 
    Z2 = @inbounds fields.Z2[i, j, k]
    D3 = @inbounds fields.D3[i, j, k] 
    B1 = @inbounds fields.B1[i, j, k]
    B2 = @inbounds fields.B2[i, j, k]

    if sum(B1+B2) > 0
        return 0.3*(bacteria_mortality(mlinB, mqB, B1) + bacteria_mortality(mlinB, mqB, B2) + phytoplankton_mortality(mlinP, mqP, P1) + phytoplankton_mortality(mlinP, mqP, P2) + zooplankton_mortality(mlinZ, mqZ, Z1) + zooplankton_mortality(mlinZ, mqZ, Z2)) - bacteria_production(μᵇ*2/2, kᴰ, y, D3, B2) / y 
    elseif sum(B1+B2) == 0
        return 0.3*(phytoplankton_mortality(mlinP, mqP, P1) + phytoplankton_mortality(mlinP, mqP, P2) + zooplankton_mortality(mlinZ, mqZ, Z1) + zooplankton_mortality(mlinZ, mqZ, Z2)) - detritus_remineralization(r/2, D3)
    end 
end

@inline function (bgc::MultiNPZBD)(i, j, k, grid, ::Val{:D4}, clock, fields)
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

    P1 = @inbounds fields.P1[i, j, k] 
    P2 = @inbounds fields.P2[i, j, k] 
    Z1 = @inbounds fields.Z1[i, j, k] 
    Z2 = @inbounds fields.Z2[i, j, k]
    D4 = @inbounds fields.D4[i, j, k] 
    B1 = @inbounds fields.B1[i, j, k]
    B2 = @inbounds fields.B2[i, j, k]

    if sum(B1+B2) > 0
        return 0.2*(bacteria_mortality(mlinB, mqB, B1) + bacteria_mortality(mlinB, mqB, B2) + phytoplankton_mortality(mlinP, mqP, P1) + phytoplankton_mortality(mlinP, mqP, P2) + zooplankton_mortality(mlinZ, mqZ, Z1) + zooplankton_mortality(mlinZ, mqZ, Z2)) - bacteria_production(μᵇ, kᴰ, y, D4, B1) / y 
    elseif sum(B1+B2) == 0
        return 0.2*(phytoplankton_mortality(mlinP, mqP, P1) + phytoplankton_mortality(mlinP, mqP, P2) + zooplankton_mortality(mlinZ, mqZ, Z1) + zooplankton_mortality(mlinZ, mqZ, Z2)) - detritus_remineralization(r, D4)
    end 
end

@inline function (bgc::MultiNPZBD)(i, j, k, grid, ::Val{:D5}, clock, fields)
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

    P1 = @inbounds fields.P1[i, j, k] 
    P2 = @inbounds fields.P2[i, j, k] 
    Z1 = @inbounds fields.Z1[i, j, k] 
    Z2 = @inbounds fields.Z2[i, j, k]
    D5 = @inbounds fields.D5[i, j, k] 
    B1 = @inbounds fields.B1[i, j, k]
    B2 = @inbounds fields.B2[i, j, k]

    if sum(B1+B2) > 0
        return 0.1*(bacteria_mortality(mlinB, mqB, B1) + bacteria_mortality(mlinB, mqB, B2) + phytoplankton_mortality(mlinP, mqP, P1) + phytoplankton_mortality(mlinP, mqP, P2) + zooplankton_mortality(mlinZ, mqZ, Z1) + zooplankton_mortality(mlinZ, mqZ, Z2)) - bacteria_production(μᵇ*2, kᴰ, y, D5, B2) / y 
    elseif sum(B1+B2) == 0
        return 0.1*(phytoplankton_mortality(mlinP, mqP, P1) + phytoplankton_mortality(mlinP, mqP, P2) + zooplankton_mortality(mlinZ, mqZ, Z1) + zooplankton_mortality(mlinZ, mqZ, Z2)) - detritus_remineralization(r, D5)
    end 
end