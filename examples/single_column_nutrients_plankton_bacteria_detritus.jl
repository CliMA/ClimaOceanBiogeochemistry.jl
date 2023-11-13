# # Nutrients, plankton, bacteria, detritus
#
# This example illustrates how to use ClimaOceanBiogeochemistry's
# `NutrientsPlanktonBacteriaDetrius` model in a single column context.

using ClimaOceanBiogeochemistry: NutrientsPlanktonBacteriaDetritus

using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Printf
using GLMakie #CairoMakie

# ## A single column grid
#
# We set up a single column grid with Nz m grid spacing that's H m deep:

H = 1000
z = (-H, 0)
Nz = 100

grid = RectilinearGrid(size = Nz; z, topology = (Flat, Flat, Bounded))

# A prescribed vertical tracer diffusivity
# 
# We define a tracer diffusivity that mixes a lot near the surface
# (in the top 50 m), and less down below.

@inline κ(x, y, z, t) = 1e-4 + 1e-2 * exp(z / 25) + 1e-2 * exp(-(z + 1000) / 50)
vertical_diffusion = VerticalScalarDiffusivity(; κ)

# We put the pieces together.
# The important line here is `biogeochemistry = NutrientsPlanktonBacteriaDetritus()`.
# We use all default parameters.

model = HydrostaticFreeSurfaceModel(; grid,
                                    velocities = PrescribedVelocityFields(),
                                    biogeochemistry = NutrientsPlanktonBacteriaDetritus(),
                                    tracers = (:N, :P, :Z, :B, :D),
                                    tracer_advection = WENO(),
                                    buoyancy = nothing,
                                    closure = vertical_diffusion)
                                    #closure = CATKEVerticalDiffusivity(),
                                    #tracers = (:b, :e),
                                    #buoyancy = BuoyancyTracer(),
                                    #boundary_conditions = (; b=b_bcs))

# ## Initial conditions
#
# We initialize the model with reasonable nutrients, detritus, and a nutrient
# mixed layer.

set!(model, P=1e-1, Z=1e-1, B=1e-2, D=5e-1, N=2)

simulation = Simulation(model, Δt=30minutes, stop_time=100days)

function progress(sim)
    @printf("Iteration: %d, time: %s, max(P): %.2e, max(Z): %.2e, max(N): %.2e, max(B): %.2e, max(D): %.2e \n",
            iteration(sim), prettytime(sim),
            maximum(model.tracers.P),
            maximum(model.tracers.Z),
            maximum(model.tracers.N),
            maximum(model.tracers.B),
            maximum(model.tracers.D))
    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

N = model.tracers.N
P = model.tracers.P
Z = model.tracers.Z
B = model.tracers.B
D = model.tracers.D

z = znodes(N)

fig = Figure()

axN = Axis(fig[1, 1], xlabel="Nutrient concentration (N)", ylabel="z (m)")
axP = Axis(fig[1, 2], xlabel="Phytoplankton concentration (P)", ylabel="z (m)")
axZ = Axis(fig[1, 3], xlabel="Zooplankton concentration (Z)", ylabel="z (m)")
axB = Axis(fig[1, 4], xlabel="Bacteria concentration (B)", ylabel="z (m)")
axD = Axis(fig[1, 5], xlabel="Detritus concentration (D)", ylabel="z (m)")

lines!(axN, interior(N, 1, 1, :), z)
lines!(axP, interior(P, 1, 1, :), z)
lines!(axZ, interior(Z, 1, 1, :), z)
lines!(axB, interior(B, 1, 1, :), z)
lines!(axD, interior(D, 1, 1, :), z)

display(fig)

filename = "nutrients_plankton_bacteria_detritus.jld2"

simulation.output_writers[:fields] = JLD2OutputWriter(model, model.tracers;
                                                      filename,
                                                      schedule = TimeInterval(1day),
                                                      overwrite_existing = true)

run!(simulation)

# ## Visualization
#
# All that's left is to visualize the results.

Pt = FieldTimeSeries(filename, "P")
Zt = FieldTimeSeries(filename, "Z")
Bt = FieldTimeSeries(filename, "B")
Dt = FieldTimeSeries(filename, "D")
Nt = FieldTimeSeries(filename, "N")

t = Pt.times
nt = length(t)
z = znodes(Pt)

fig = Figure(resolution=(1200, 600))

axN = Axis(fig[1, 1], ylabel="z (m)", xlabel="Nutrient concentration (mmol)")
axP = Axis(fig[1, 2], ylabel="z (m)", xlabel="Phytoplankton concentration (mmol)")
axZ = Axis(fig[1, 3], ylabel="z (m)", xlabel="Zooplankton concentration (mmol)")
axB = Axis(fig[1, 4], ylabel="z (m)", xlabel="Bacteria concentration (mmol)")
axD = Axis(fig[1, 5], ylabel="z (m)", xlabel="Detritus concentration (mmol)")

xlims!(axP, -0.1, 1)
xlims!(axZ, -0.1, 0.5)
xlims!(axB, -0.1, 0.5)
xlims!(axD, -0.1, 0.5)

slider = Slider(fig[2, 1:5], range=1:nt, startvalue=1)
n = slider.value

title = @lift @sprintf("Equilibrium biogeochemistry at t = %d days", t[$n] / day)
Label(fig[0, 1:5], title)

Nn = @lift interior(Nt[$n], 1, 1, :)
Pn = @lift interior(Pt[$n], 1, 1, :)
Zn = @lift interior(Zt[$n], 1, 1, :)
Bn = @lift interior(Bt[$n], 1, 1, :)
Dn = @lift interior(Dt[$n], 1, 1, :)

lines!(axP, Pn, z)
lines!(axZ, Zn, z)
lines!(axD, Dn, z)
lines!(axB, Bn, z)
lines!(axN, Nn, z)

record(fig, "nutrients_plankton_bacteria_detritus.mp4", 1:nt, framerate=24) do nn
    n[] = nn
end
nothing #hide

# ![](nutrients_plankton_bacteria_detritus.mp4)
