# # test POP remineralization formulations in 2D settings
#
# The BGC tracers in our model are advected, diffuse, precipitate, and dissolve 

# We use Oceananigans'`Forcing` abstraction to implement the P cycle dynamics
#
# This example demonstrates
#   * How to set time-dependent boundary conditions.
#   * How to use the `TimeStepWizard` to adapt the simulation time-step.
#   * How to use `Average` to diagnose spatial averages of model fields.

using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients
using GLMakie
using Printf
using Statistics

using Oceananigans
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Oceananigans.Units

#################################### Grid ####################################
# Parameters
Nx = 100
Nz = 100
Lx = 10kilometers   # m
Lz = 1000            # m

# We use a two-dimensional grid, with a `Flat` `y`-direction:
grid = RectilinearGrid(size = (Nx, Nz),
                       x = (0, Lx),
                       z = (-Lz, 0),
                       topology=(Bounded, Flat, Bounded))

#################################### Boundary conditions ####################################

# Compute piston velocity `wp` from restoring time-scale τ★
Δz = zspacings(grid, Center())
τ★ = 1days  
wp = Δz / τ★ # piston velocity = 10 m/day

# Compute buoyancy flux Jᵇ
M² = 1e-12         # s⁻², squared buoyancy frequency
@inline function Jᵇ(x, t, b, p)
    b★ = p.M² * x   # m s⁻², reference buoyancy
    return p.wp * (b★ - b) # m² s⁻³
end
top_buoyancy_bc = FluxBoundaryCondition(Jᵇ, field_dependencies=:b, parameters=(; M², wp))
b_bcs = FieldBoundaryConditions(top=top_buoyancy_bc)

# Wind stress boundary condition
τ₀ = 1e-4           # m² s⁻²
@inline step(x, xc, dx) = (1 + tanh((x - xc) / dx)) / 2 # between 0 and 1
@inline τx(x, t, p) = - p.τ₀ * step(x, p.Lx/2, p.Lx/5)
x_wind_stress = FluxBoundaryCondition(τx, parameters=(; τ₀, Lx))
u_bcs = FieldBoundaryConditions(top=x_wind_stress)

#################################### Forcing ####################################
# buoyancy increase with depth: downwelling
buoyancy_gradient = 1e-6 # m s⁻² m⁻¹ # positive: increase with depth
surface_buoyancy = -1e-6    # m s⁻² (at z=0)

target_buoyancy = LinearTarget{:z}(intercept=surface_buoyancy, gradient=buoyancy_gradient)
Δx = Lx/Nx
x_mask = GaussianMask{:x}(center=Δx/2, width=Δx/2)
b_sponge = Relaxation(rate=1/1days, mask=x_mask, target=target_buoyancy)

#################################### Model ####################################

background_diffusivity = ScalarDiffusivity(ν=1e-4, κ=1e-4)
catke = CATKEVerticalDiffusivity()

model = HydrostaticFreeSurfaceModel(; grid,
                                    buoyancy = BuoyancyTracer(),
                                    biogeochemistry = CarbonAlkalinityNutrients(; grid,
                                                                        maximum_net_community_production_rate  = 5e-2/day,
                                                                        fraction_of_particulate_export = 1),
                                    coriolis = FPlane(; f=1e-4),
                                    closure = (catke, background_diffusivity),
                                    tracers = (:b, :e, :DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe),
                                    # tracers = (:b, :e), 
                                    # forcing=(;b=b_sponge),
                                    boundary_conditions = (; u=u_bcs)) #, b=b_bcs))

N² = 1e-5           # s⁻²
bᵢ(x, z) = N² * z + M² * x

# Biogeochemical
N₀ = 1e-2 # Surface nutrient concentration
D₀ = 1e-3 # Surface detritus concentration

Nᵢ(x,z) = N₀ * (-z / 115).^0.84
Dᵢ(x,z) = D₀ * (-z / 115).^-0.84

set!(model, b=bᵢ, DIC=2.1, ALK=2.35, NO₃=1e-2, PO₄=Nᵢ, DOP=Dᵢ, POP=Dᵢ, Fe = 1e-6) # mol PO₄ m⁻³
# set!(model, b=bᵢ) 

# TODO: When Δt is large (>30 seconds), the model reaches NaN quickly
simulation = Simulation(model; Δt = 20seconds, stop_time=100days)

# progress(sim) = @printf("Iteration: %d, time: %s, total(P): %.2e\n",
#                         iteration(sim), prettytime(sim), 
#                         sum(model.tracers.PO₄)+sum(model.tracers.DOP)+sum(model.tracers.POP))
progress(sim) = @printf("Iteration: %d, time: %s\n",
                        iteration(sim), prettytime(sim))

add_callback!(simulation, progress, IterationInterval(100))

outputs = (u = model.velocities.u,
            w = model.velocities.w,
            b = model.tracers.b,
            PO₄= model.tracers.PO₄,
            DOP = model.tracers.DOP,
            POP = model.tracers.POP,
            avg_PO₄ = Average(model.tracers.PO₄, dims=(1, 2)),
            avg_DOP = Average(model.tracers.DOP, dims=(1, 2)),
            avg_POP = Average(model.tracers.POP, dims=(1, 2)))

simulation.output_writers[:simple_output] =
    JLD2OutputWriter(model, outputs,
                     schedule = TimeInterval(30minutes),
                     filename = "Remin_2D",
                     overwrite_existing = true)

run!(simulation)

############################ Visualizing the solution ############################

filepath = simulation.output_writers[:simple_output].filepath

u_timeseries = FieldTimeSeries(filepath, "u")
w_timeseries = FieldTimeSeries(filepath, "w")
times = w_timeseries.times
xw, yw, zw = nodes(w_timeseries)
b_timeseries = FieldTimeSeries(filepath, "b")
xi, yi, zi = nodes(b_timeseries)

PO4_timeseries = FieldTimeSeries(filepath, "PO₄")
avg_PO4_timeseries = FieldTimeSeries(filepath, "avg_PO₄")
POP_timeseries = FieldTimeSeries(filepath, "POP")
avg_POP_timeseries = FieldTimeSeries(filepath, "avg_POP")
# xi, yi, zi = nodes(avg_PO4_timeseries)

n = Observable(1)
title = @lift @sprintf("t = %d days", times[$n] / day)

uₙ = @lift interior(u_timeseries[$n], :, 1, :)
wₙ = @lift interior(w_timeseries[$n], :, 1, :)
bₙ = @lift interior(b_timeseries[$n], :, 1, :)

PO4ₙ = @lift interior(PO4_timeseries[$n], :, 1, :)
avg_PO4ₙ = @lift interior(avg_PO4_timeseries[$n], 1, 1, :)

POPₙ = @lift interior(POP_timeseries[$n], :, 1, :)
avg_POPₙ = @lift interior(avg_POP_timeseries[$n], 1, 1, :)

fig = Figure(size = (1200, 1200))

ax_u = Axis(fig[2, 1]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_u = heatmap!(ax_u, xw, zw, uₙ; colorrange = (-0.02, 0.02), colormap = :balance)
Colorbar(fig[2, 2], hm_u; label = "u", flipaxis = false)

ax_w = Axis(fig[2, 3]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_w = heatmap!(ax_w, xw, zw, wₙ; colorrange = (-0.0002, 0.0002), colormap = :balance)
Colorbar(fig[2, 4], hm_w; label = "w", flipaxis = false)

ax_b = Axis(fig[2, 5]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_b = heatmap!(ax_b, xi, zi, bₙ; colorrange = (-0.005, 0), colormap = :matter)
Colorbar(fig[2, 6], hm_b; label = "b", flipaxis = false)

ax_PO4 = Axis(fig[3, 1]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_PO4 = heatmap!(ax_PO4, xw, zw, PO4ₙ; colormap = :matter, colorrange = (0.0, 0.05))
Colorbar(fig[3, 2], hm_PO4; label = "[PO₄]", flipaxis = false)

ax_POP = Axis(fig[3, 3]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_POP = heatmap!(ax_POP, xw, zw, POPₙ; colormap = :matter, colorrange = (0.0, 0.02))
Colorbar(fig[3, 4], hm_POP; label = "[POP]", flipaxis = false)

ax_avg_PO4 = Axis(fig[4, 1]; xlabel = "[PO₄] (μM)", ylabel = "z (m)", yaxisposition = :right)
xlims!(ax_avg_PO4, 0.0, 0.08)
lines!(ax_avg_PO4, avg_PO4ₙ, zi)

ax_avg_POP = Axis(fig[4, 3]; xlabel = "[POP] (μM)", ylabel = "z (m)", yaxisposition = :right)
xlims!(ax_avg_POP, 0.0, 0.01)
lines!(ax_avg_POP, avg_POPₙ, zi)

fig[1, 1:6] = Label(fig, title, tellwidth=false)

# And, finally, we record a movie.
frames = 1:length(times)
record(fig, "Remin_2D.mp4", frames, framerate=24) do i
    n[] = i
end
nothing #hide
