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

# We use a two-dimensional grid with 64² points, 3² halo points for high-order advection,
# 2 m grid spacing, and a `Flat` `y`-direction:

grid = RectilinearGrid(size=(64, 64), extent=(512, 512), halo=(3, 3), topology=(Periodic, Flat, Bounded))

# ## Boundary conditions
#
# We impose a surface buoyancy flux that's initially constant and then decays to zero,

buoyancy_flux(x, t, params) = params.initial_buoyancy_flux * exp(-t^4 / (24 * params.shut_off_time^4))

buoyancy_flux_parameters = (initial_buoyancy_flux = 1e-8, # m² s⁻³
                                    shut_off_time = 2hours)

buoyancy_flux_bc = FluxBoundaryCondition(buoyancy_flux, parameters = buoyancy_flux_parameters)

# The buoyancy flux effectively shuts off after 6 hours of simulation time.
# The initial condition and bottom boundary condition impose the constant buoyancy gradient

N² = 1e-4 # s⁻²
buoyancy_gradient_bc = GradientBoundaryCondition(N²)

# In summary, the buoyancy boundary conditions impose a destabilizing flux
# at the top and a stable buoyancy gradient at the bottom:

buoyancy_bcs = FieldBoundaryConditions(top = buoyancy_flux_bc, bottom = buoyancy_gradient_bc)

# ## Phosphorus model with parameters

biogeochemistry = CarbonAlkalinityNutrients(; grid,
                                            maximum_net_community_production_rate  = 5e-3/day,
                                            dissolved_organic_phosphorus_remin_timescale  = 0.05/day,
                                            particulate_organic_phosphorus_sinking_velocity = -1.0/day)
                            
# We tell `Forcing` that our CAN model depends
# on the plankton concentration `P` and the chosen parameters,

# P_dynamics = Forcing(growing_and_grazing, field_dependencies = :P,
#                             parameters = plankton_dynamics_parameters)

# ## The model
#
# The name "`P`" for phytoplankton is specified in the
# constructor for `NonhydrostaticModel`. We additionally specify a fifth-order
# advection scheme, third-order Runge-Kutta time-stepping, isotropic viscosity and diffusivities,
# and Coriolis forces appropriate for planktonic convection at mid-latitudes on Earth.

################## TODO: How to set forcing? ################
model = NonhydrostaticModel(; grid,
                            advection = UpwindBiasedFifthOrder(),
                            biogeochemistry,
                            tracers = (:b, :DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe),
                            timestepper = :RungeKutta3,
                            closure = ScalarDiffusivity(ν=1e-4, κ=1e-4),
                            coriolis = FPlane(f=1e-4),
                            buoyancy = BuoyancyTracer(),
                            # forcing = (; P=P_dynamics),
                            boundary_conditions = (; b=buoyancy_bcs))

# ## Initial condition
#
# For buoyancy, we use a stratification that's mixed near the surface and
# linearly stratified below, superposed with surface-concentrated random noise.

mixed_layer_depth = 32 # m
stratification(z) = z < -mixed_layer_depth ? N² * z : - N² * mixed_layer_depth
noise(z) = 1e-4 * N² * grid.Lz * randn() * exp(z / 4)
initial_buoyancy(x, z) = stratification(z) + noise(z)

set!(model, b=initial_buoyancy, DIC=2.1, ALK=2.35, NO₃=3.6e-2, PO₄=2e-3, DOP=1e-3, POP=1e-3, Fe = 1e-6)

# ## Simulation with adaptive time-stepping, logging, and output
#
# We build a simulation

simulation = Simulation(model, Δt=2minutes, stop_time=20days)

# with a `TimeStepWizard` that limits the
# time-step to 2 minutes, and adapts the time-step such that CFL
# (Courant-Freidrichs-Lewy) number hovers around `1.0`,

conjure_time_step_wizard!(simulation, cfl=1.0, max_Δt=2minutes)

# We also add a callback that prints the progress of the simulation,

progress(sim) = @printf("Iteration: %d, time: %s, total(P): %.2e\n",
                        iteration(sim), prettytime(sim), 
                        sum(model.tracers.PO₄)+sum(model.tracers.DOP)+sum(model.tracers.POP))


add_callback!(simulation, progress, IterationInterval(100))

# and a basic `JLD2OutputWriter` that writes velocities and both
# the two-dimensional and horizontally-averaged concentration

outputs = (w = model.velocities.w,
            PO₄= model.tracers.PO₄,
            DOP = model.tracers.DOP,
            POP = model.tracers.POP,
            avg_PO₄ = Average(model.tracers.PO₄, dims=(1, 2)),
            avg_DOP = Average(model.tracers.DOP, dims=(1, 2)),
            avg_POP = Average(model.tracers.POP, dims=(1, 2)))

simulation.output_writers[:simple_output] =
    JLD2OutputWriter(model, outputs,
                     schedule = TimeInterval(30minutes),
                     filename = "convecting_plankton.jld2",
                     overwrite_existing = true)

run!(simulation)

# Notice how the time-step is reduced at early times, when turbulence is strong,
# and increases again towards the end of the simulation when turbulence fades.

# ## Visualizing the solution


filepath = simulation.output_writers[:simple_output].filepath

w_timeseries = FieldTimeSeries(filepath, "w")
PO4_timeseries = FieldTimeSeries(filepath, "PO₄")
avg_PO4_timeseries = FieldTimeSeries(filepath, "avg_PO₄")
DOP_timeseries = FieldTimeSeries(filepath, "DOP")
avg_DOP_timeseries = FieldTimeSeries(filepath, "avg_DOP")
POP_timeseries = FieldTimeSeries(filepath, "POP")
avg_POP_timeseries = FieldTimeSeries(filepath, "avg_POP")

times = w_timeseries.times
# buoyancy_flux_time_series = [buoyancy_flux(0, t, buoyancy_flux_parameters) for t in times]
# nothing #hide

# and then we construct the ``x, z`` grid,

# xw, yw, zw = nodes(w_timeseries)
xi, yi, zi = nodes(PO4_timeseries) # inorganic
xd, yd, zd = nodes(DOP_timeseries) # dissolved
xp, yp, zp = nodes(POP_timeseries) # particulate
nothing #hide

# Finally, we animate 

@info "Making a movie..."

n = Observable(1)
title = @lift @sprintf("t = %s", prettytime(times[$n]))

# wₙ = @lift interior(w_timeseries[$n], :, 1, :)
PO4ₙ = @lift interior(PO4_timeseries[$n], :, 1, :)
DOPₙ = @lift interior(DOP_timeseries[$n], :, 1, :)
POPₙ = @lift interior(POP_timeseries[$n], :, 1, :)
avg_PO4ₙ = @lift interior(avg_PO4_timeseries[$n], 1, 1, :)
avg_DOPₙ = @lift interior(avg_DOP_timeseries[$n], 1, 1, :)
avg_POPₙ = @lift interior(avg_POP_timeseries[$n], 1, 1, :)

PO4_lims = (0.0, 0.005)
DOP_lims = (0.0, 0.005)
POP_lims = (0.0, 0.005)

fig = Figure(size = (1000, 1000))

ax_PO4 = Axis(fig[2, 2]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
ax_DOP = Axis(fig[3, 2]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
ax_POP = Axis(fig[4, 2]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
ax_avg_PO4 = Axis(fig[2, 3]; xlabel = "[PO₄] (μM)", ylabel = "z (m)", yaxisposition = :right)
ax_avg_DOP = Axis(fig[3, 3]; xlabel = "[DOP] (μM)", ylabel = "z (m)", yaxisposition = :right)
ax_avg_POP = Axis(fig[4, 3]; xlabel = "[POP] (μM)", ylabel = "z (m)", yaxisposition = :right)
xlims!(ax_avg_PO4, 0.0, 0.005)
xlims!(ax_avg_DOP, 0.0, 0.005)
xlims!(ax_avg_POP, 0.0, 0.005)

fig[1, 1:3] = Label(fig, title, tellwidth=false)

hm_PO4 = heatmap!(ax_PO4, xi, zi, PO4ₙ; colormap = :matter, colorrange = PO4_lims)
Colorbar(fig[2, 1], hm_PO4; label = "[PO₄]", flipaxis = false)
hm_DOP = heatmap!(ax_DOP, xd, zd, DOPₙ; colormap = :matter, colorrange = DOP_lims)
Colorbar(fig[3, 1], hm_DOP; label = "[DOP]", flipaxis = false)
hm_POP = heatmap!(ax_POP, xp, zp, POPₙ; colormap = :matter, colorrange = POP_lims)
Colorbar(fig[4, 1], hm_POP; label = "[POP]", flipaxis = false)

lines!(ax_avg_PO4, avg_PO4ₙ, zi)
lines!(ax_avg_DOP, avg_DOPₙ, zd)
lines!(ax_avg_POP, avg_POPₙ, zp)

current_figure() #hide
fig

# And, finally, we record a movie.

frames = 1:length(times)

record(fig, "convecting_plankton.mp4", frames, framerate=8) do i
    n[] = i
end
nothing #hide

# ![](convecting_plankton.mp4)
