using GLMakie
using Oceananigans
using Oceananigans.Units
using Printf

using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity

# Parameters
Nz = 50
Lz = 400
f = 1e-4   # Coriolis parameter
N² = 1e-5  
Qᵇ = +1e-8 # cooling buoyancy flux
Qᵘ = -2e-4 # drives flow in +x direction

grid = RectilinearGrid(size=Nz, z=(-Lz, 0), topology=(Flat, Flat, Bounded))

@info "Using $grid"

b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵇ))
u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵘ))

vertical_mixing = CATKEVerticalDiffusivity()

model = HydrostaticFreeSurfaceModel(; grid,
                                    coriolis = FPlane(; f),
                                    closure = vertical_mixing,
                                    tracers = (:b, :e),
                                    buoyancy = BuoyancyTracer(),
                                    boundary_conditions = (; b=b_bcs, u=u_bcs))
                                    
@info "Build $model"

bᵢ(x, y, z) = N² * z
set!(model, b=bᵢ, e=1e-6)

simulation = Simulation(model, Δt=20minute, stop_time=8days)

simulation.output_writers[:fields] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers),
                     schedule = TimeInterval(1hour),
                     filename = "windy_convection.jld2",
                     overwrite_existing = true)

progress(sim) = @info string("Iter: ", iteration(sim), " t: ", prettytime(sim))
simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

@info "Running a simulation of $model..."

run!(simulation)

#####
##### Visualize
#####

filepath = "windy_convection.jld2"

bt = FieldTimeSeries(filepath, "b")
ut = FieldTimeSeries(filepath, "u")
vt = FieldTimeSeries(filepath, "v")
et = FieldTimeSeries(filepath, "e")

z = znodes(bt)
times = bt.times
Nt = length(times)
ulim = 1.1 * maximum(abs, ut)
elim = 1.1 * maximum(abs, et[Nt])

fig = Figure(resolution=(1200, 800))

slider = Slider(fig[2, 1:2], range=1:Nt, startvalue=1)
n = slider.value

buoyancy_label = @lift "Buoyancy at t = " * prettytime(times[$n])
velocities_label = @lift "Velocities at t = " * prettytime(times[$n])
TKE_label = @lift "Turbulent kinetic energy t = " * prettytime(times[$n])

ax_b = Axis(fig[1, 1], xlabel=buoyancy_label, ylabel="z (m)")
ax_u = Axis(fig[1, 2], xlabel=velocities_label, ylabel="z (m)")
ax_e = Axis(fig[1, 3], xlabel=TKE_label, ylabel="z (m)")

xlims!(ax_b, -grid.Lz * N², 0)
xlims!(ax_u, -ulim, ulim)
xlims!(ax_e, 0, elim)

colors = [:black, :blue, :red, :orange]

bn = @lift interior(bt[$n], 1, 1, :)
un = @lift interior(ut[$n], 1, 1, :)
vn = @lift interior(vt[$n], 1, 1, :)
en = @lift interior(et[$n], 1, 1, :)

lines!(ax_b, bn, z)
lines!(ax_u, un, z, label="u")
lines!(ax_u, vn, z, label="v")
lines!(ax_e, en, z)

axislegend(ax_u, position=:rb)

display(fig)

record(fig, "windy_convection.mp4", 1:Nt, framerate=24) do nn
    n[] = nn
end

