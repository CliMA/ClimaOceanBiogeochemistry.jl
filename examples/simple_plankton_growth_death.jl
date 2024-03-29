using Oceananigans
using Oceananigans.Units
using Oceananigans
using Oceananigans.Grids: znode
using Oceananigans.Biogeochemistry: AbstractBiogeochemistry
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity

using Printf

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers

"""
    SimplePlanktonGrowthDeath(; growth_rate = 1/day
                                shortwave_attenuation_scale = 20.0
                                mortality_rate = 0.1/day)

Return a single-tracer biogeochemistry model for growing and dying plankton.
"""
Base.@kwdef struct SimplePlanktonGrowthDeath{FT} <: AbstractBiogeochemistry
    growth_rate :: FT = 1/day
    shortwave_attenuation_scale :: FT = 20.0
    mortality_rate :: FT = 0.1/day
end

@inline required_biogeochemical_tracers(::SimplePlanktonGrowthDeath) = tuple(:P)

const c = Center()

@inline function (bgc::SimplePlanktonGrowthDeath)(i, j, k, grid, ::Val{:P}, clock, fields)
   μ₀ = bgc.growth_rate
   λ = bgc.shortwave_attenuation_scale
   m = bgc.mortality_rate
   P = @inbounds fields.P[i, j, k]
   z = znode(i, j, k, grid, c, c, c)
   return (μ₀ * exp(z / λ) - m) * P
end

# We set up the model

grid = RectilinearGrid(size = 64,
                       z = (-256meters, 0),
                       topology = (Flat, Flat, Bounded))

Qᵇ(t) = ifelse(t < 4days, 1e-7, 0.0)
b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵇ))

model = HydrostaticFreeSurfaceModel(; grid,
                                    biogeochemistry = SimplePlanktonGrowthDeath(),
                                    closure = CATKEVerticalDiffusivity(),
                                    tracers = (:b, :e),
                                    buoyancy = BuoyancyTracer(),
                                    boundary_conditions = (; b=b_bcs))

N² = 1e-5 # s⁻²
bᵢ(z) = N² * z
set!(model, b=bᵢ, P=1e-2, e=1e-6)

simulation = Simulation(model, Δt=10minutes, stop_time=8days)

progress(sim) = @printf("Iteration: %d, time: %s, max(P): %.2e \n",
                        iteration(sim), prettytime(sim), maximum(model.tracers.P))

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

outputs = merge(model.velocities, model.tracers)
filename = "simple_plankton_growth_death.jld2"

simulation.output_writers[:fields] = JLD2OutputWriter(model, outputs;
                                                      filename,
                                                      schedule = TimeInterval(20minutes),
                                                      overwrite_existing = true)

run!(simulation)

# Now we load the saved output and plot

using CairoMakie

bt = FieldTimeSeries(filename, "b")
et = FieldTimeSeries(filename, "e")
Pt = FieldTimeSeries(filename, "P")

t = bt.times
Nt = length(t)
z = znodes(bt)

fig = Figure(size=(800, 400))

axb = Axis(fig[1, 1], ylabel="z (m)", xlabel="Buoyancy (m² s⁻³)")
axe = Axis(fig[1, 2], ylabel="z (m)", xlabel="Turbulent kinetic energy (m² s²)")
axP = Axis(fig[1, 3], ylabel="z (m)", xlabel="Plankton concentration")

xlims!(axe, -1e-5, 1e-3)
xlims!(axP, 0, 0.1)

n = Observable(1)

title = @lift @sprintf("Convecting plankton at t = %d days", t[$n] / day)
Label(fig[0, 1:3], title)

bn = @lift interior(bt[$n], 1, 1, :)
en = @lift interior(et[$n], 1, 1, :)
Pn = @lift interior(Pt[$n], 1, 1, :)

lines!(axb, bn, z)
lines!(axe, en, z)
lines!(axP, Pn, z)

fig

record(fig, "simple_plankton_growth_death.mp4", 1:Nt, framerate=24) do nn
    n[] = nn
end

# ![](simple_plankton_growth_death.mp4)
