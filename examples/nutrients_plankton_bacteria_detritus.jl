using Oceananigans
using Oceananigans.Units
using Oceananigans
using Oceananigans.Grids: znode
using Oceananigans.Biogeochemistry: AbstractBiogeochemistry
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity

using Printf

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers

"""
    NutrientsPlanktonBacteriaDetritus(; kw...)

Return a four-tracer biogeochemistry model for growing and dying plankton.
"""
Base.@kwdef struct NutrientsPlanktonBacteriaDetritus{FT} <: AbstractBiogeochemistry
    maximum_plankton_growth_rate :: FT = 1/day
    maximum_bacteria_growth_rate :: FT = 1/day
    bacteria_yield :: FT               = 0.2
    quadratic_mortality_rate :: FT     = 1/day # m³/mmol/day
    nutrient_half_saturation :: FT     = 0.1   # mmol m⁻³
    detritus_half_saturation :: FT     = 0.1   # mmol m⁻³
    PAR_half_saturation :: FT          = 10.0  # W m⁻²
    PAR_attenuation_scale :: FT        = 25.0  # m
    detritus_sinking_speed :: FT       = 10/day # m s⁻¹
end

@inline required_biogeochemical_tracers(::NutrientsPlanktonBacteriaDetritus) = (:N, :P, :B, :D)

const c = Center()

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

    return bacteria_mortality(m, B) + plankton_mortality(m, B) - bacteria_production(μᵇ, kᴰ, D, B) / y
end

#####
##### Set up the model
#####

grid = RectilinearGrid(size = 64,
                       z = (-256, 0),
                       topology = (Flat, Flat, Bounded))

Qᵇ(x, y, t) = ifelse(t < 4days, 1e-7, 0.0)
b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵇ))

model = HydrostaticFreeSurfaceModel(; grid,
                                    biogeochemistry = NutrientsPlanktonBacteriaDetritus(),
                                    closure = CATKEVerticalDiffusivity(),
                                    tracers = (:b, :e),
                                    buoyancy = BuoyancyTracer(),
                                    boundary_conditions = (; b=b_bcs))

N² = 1e-5 # s⁻²
bᵢ(x, y, z) = N² * z

set!(model, b=bᵢ, P=1e-1, B=1e-1, D=1e-1, N=1e1, e=1e-6)

simulation = Simulation(model, Δt=10minutes, stop_iteration=1200)

progress(sim) = @printf("Iteration: %d, time: %s, max(P): %.2e, max(N): %.2e, max(B): %.2e, max(D): %.2e \n",
                        iteration(sim), prettytime(sim),
                        maximum(model.tracers.P),
                        maximum(model.tracers.N),
                        maximum(model.tracers.B),
                        maximum(model.tracers.D))

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

outputs = merge(model.velocities, model.tracers)
filename = "simple_plankton_growth_death.jld2"

simulation.output_writers[:fields] = JLD2OutputWriter(model, outputs;
                                                      filename,
                                                      schedule = IterationInterval(10),
                                                      overwrite_existing = true)

run!(simulation)

using GLMakie

bt = FieldTimeSeries(filename, "b")
et = FieldTimeSeries(filename, "e")
Pt = FieldTimeSeries(filename, "P")
Bt = FieldTimeSeries(filename, "B")
Dt = FieldTimeSeries(filename, "D")
Nt = FieldTimeSeries(filename, "N")

t = bt.times
nt = length(t)
z = znodes(bt)

fig = Figure(resolution=(1200, 600))

axb = Axis(fig[1, 1], ylabel="z (m)", xlabel="Buoyancy (m² s⁻³)")
axe = Axis(fig[1, 2], ylabel="z (m)", xlabel="Turbulent kinetic energy (m² s²)")
axP = Axis(fig[1, 3], ylabel="z (m)", xlabel="Concentration")
axN = Axis(fig[1, 4], ylabel="z (m)", xlabel="Nutrient concentration")

xlims!(axe, -1e-5, 1e-3)
xlims!(axP, 0, 0.2)

slider = Slider(fig[2, 1:4], range=1:nt, startvalue=1)
n = slider.value

title = @lift @sprintf("Convecting plankton at t = %d hours", t[$n] / hour)
Label(fig[0, 1:4], title)

bn = @lift interior(bt[$n], 1, 1, :)
en = @lift interior(et[$n], 1, 1, :)
Pn = @lift interior(Pt[$n], 1, 1, :)
Bn = @lift interior(Bt[$n], 1, 1, :)
Dn = @lift interior(Dt[$n], 1, 1, :)
Nn = @lift interior(Nt[$n], 1, 1, :)

lines!(axb, bn, z)
lines!(axe, en, z)
lines!(axP, Pn, z, label="Plankton")
lines!(axP, Dn, z, label="Detritus")
lines!(axP, Bn, z, label="Bacteria")
axislegend(axP)
lines!(axN, Nn, z)

display(fig)

# record(fig, "simple_plankton_growth_death.mp4", 1:Nt, framerate=24) do nn
#     n[] = nn
# end

