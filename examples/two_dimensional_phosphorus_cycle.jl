# # test POP remineralization formulations in 2D settings
#
# The BGC tracers in our model are advected, diffuse, precipitate, and dissolve 

# We use Oceananigans'`Forcing` abstraction to implement the P cycle dynamics
#
# This example demonstrates
#   * How to set time-dependent boundary conditions.
#   * How to use the `TimeStepWizard` to adapt the simulation time-step.
#   * How to use `Average` to diagnose spatial averages of model fields.

using Oceananigans
using Oceananigans.Units: minutes, hour, hours, day, days

using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients
using CairoMakie
using Printf
using Statistics

# We use a two-dimensional grid with 100² points,
# 10 m grid spacing, and a `Flat` `y`-direction:

grid = RectilinearGrid(size=(100, 100), extent=(1000, 1000), topology=(Bounded, Flat, Bounded))

#################################### Boundary conditions ####################################

velocity_bcs = FieldBoundaryConditions(top=ValueBoundaryCondition(0.0),
                                    bottom=ValueBoundaryCondition(0.0))

################################# The model #################################

biogeochemistry = CarbonAlkalinityNutrients(; grid,
                                            maximum_net_community_production_rate  = 5e-2/day,
                                            fraction_of_particulate_export  = 1,
                                            martin_curve_exponent = 2,
                                            particulate_organic_phosphorus_sinking_velocity = -1.0/day)
                            

model = NonhydrostaticModel(; grid,
                            advection = UpwindBiasedFifthOrder(),
                            biogeochemistry,
                            tracers = (:b, :DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe),
                            timestepper = :RungeKutta3,
                            closure = ScalarDiffusivity(; ν=1e-4, κ=1e-4))
                            boundary_conditions = (u=velocity_bcs, w=velocity_bcs),
                            #buoyancy = BuoyancyTracer())

################################# Initial condition #################################

# Biogeochemical
N₀ = 1e-2 # Surface nutrient concentration
D₀ = 1e-3 # Surface detritus concentration

Nᵢ(x,z) = N₀ * (-z / 115).^0.84
Dᵢ(x,z) = D₀ * (-z / 115).^-0.84

# Physical
bᵢ(x, z) = 1e-4 * z

u, v, w = model.velocities

uᵢ = rand(size(u)...)
wᵢ = rand(size(w)...)

uᵢ .-= mean(uᵢ)
wᵢ .-= mean(wᵢ)

set!(model, u=uᵢ, w=wᵢ, b=bᵢ, DIC=2.1, ALK=2.35, NO₃=1e-2, PO₄=Nᵢ, DOP=0, POP=Dᵢ, Fe = 1e-6) # mol PO₄ m⁻³

############################## We build a simulation ##############################

simulation = Simulation(model, Δt=5minutes, stop_time=10days)

# with a `TimeStepWizard` that limits the
# time-step to 2 minutes, and adapts the time-step such that CFL
# (Courant-Freidrichs-Lewy) number hovers around `1.0`,

#conjure_time_step_wizard!(simulation, cfl=1.0, max_Δt=2minutes)

# We also add a callback that prints the progress of the simulation,

# progress(sim) = @printf("Iteration: %d, time: %s, total(P): %.2e\n",
#                         iteration(sim), prettytime(sim), 
#                         sum(model.tracers.PO₄)+sum(model.tracers.DOP)+sum(model.tracers.POP))

# add_callback!(simulation, progress, IterationInterval(100))


# and a basic `JLD2OutputWriter` that writes velocities and both
# the two-dimensional and horizontally-averaged tracer concentrations

outputs = (u = model.velocities.u,
            w = model.velocities.w,
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
PO4_timeseries = FieldTimeSeries(filepath, "PO₄")
avg_PO4_timeseries = FieldTimeSeries(filepath, "avg_PO₄")
POP_timeseries = FieldTimeSeries(filepath, "POP")
avg_POP_timeseries = FieldTimeSeries(filepath, "avg_POP")

times = w_timeseries.times

# and then we construct the ``x, z`` grid,

# xw, yw, zw = nodes(w_timeseries)
xi, yi, zi = nodes(PO4_timeseries) # inorganic
# xd, yd, zd = nodes(DOP_timeseries) # dissolved
xp, yp, zp = nodes(POP_timeseries) # particulate
nothing #hide

# Finally, we animate 

@info "Making a movie..."

n = Observable(1)
title = @lift @sprintf("t = %s", prettytime(times[$n]))

uₙ = @lift interior(u_timeseries[$n], :, 1, :)
wₙ = @lift interior(w_timeseries[$n], :, 1, :)

PO4ₙ = @lift interior(PO4_timeseries[$n], :, 1, :)
avg_PO4ₙ = @lift interior(avg_PO4_timeseries[$n], 1, 1, :)

POPₙ = @lift interior(POP_timeseries[$n], :, 1, :)
avg_POPₙ = @lift interior(avg_POP_timeseries[$n], 1, 1, :)
# DOPₙ = @lift interior(DOP_timeseries[$n], :, 1, :)
# avg_DOPₙ = @lift interior(avg_DOP_timeseries[$n], 1, 1, :)


PO4_lims = (0.0, 0.01)
# DOP_lims = (0.0, 0.005)
POP_lims = (0.0, 0.01)

fig = Figure(size = (1000, 1000))

ax_u = Axis(fig[2, 2]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_u = heatmap!(ax_u, xi, zi, uₙ; colormap = :balance)
Colorbar(fig[2, 1], hm_u; label = "u", flipaxis = false)

ax_w = Axis(fig[2, 3]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_w = heatmap!(ax_w, xi, zi, wₙ; colormap = :balance)

ax_PO4 = Axis(fig[3, 2]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
ax_POP = Axis(fig[4, 2]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)

ax_avg_PO4 = Axis(fig[3, 3]; xlabel = "[PO₄] (μM)", ylabel = "z (m)", yaxisposition = :right)
ax_avg_POP = Axis(fig[4, 3]; xlabel = "[POP] (μM)", ylabel = "z (m)", yaxisposition = :right)
xlims!(ax_avg_PO4, 0.0, 0.01)
xlims!(ax_avg_POP, 0.0, 0.01)

fig[1, 1:3] = Label(fig, title, tellwidth=false)

hm_PO4 = heatmap!(ax_PO4, xi, zi, PO4ₙ; colormap = :matter, colorrange = PO4_lims)
Colorbar(fig[3, 1], hm_PO4; label = "[PO₄]", flipaxis = false)
hm_POP = heatmap!(ax_POP, xp, zp, POPₙ; colormap = :matter, colorrange = POP_lims)
Colorbar(fig[4, 1], hm_POP; label = "[POP]", flipaxis = false)

lines!(ax_avg_PO4, avg_PO4ₙ, zi)
lines!(ax_avg_POP, avg_POPₙ, zp)

current_figure() #hide
fig

# And, finally, we record a movie.
frames = 1:length(times)
record(fig, "Remin_2D.mp4", frames, framerate=8) do i
    n[] = i
end
nothing #hide
