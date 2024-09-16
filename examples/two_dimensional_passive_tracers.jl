# # test a 2D physical setting with passive tracers

using GLMakie
using Printf
using Statistics

using Oceananigans
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Oceananigans.Units

#################################### Grid ####################################
# Parameters
Nx = 100 
Nz = 200
Lx = 1000kilometers   # m
Lz = 2000           # m

# We use a two-dimensional grid, with a `Flat` `y`-direction:
grid = RectilinearGrid(size = (Nx, Nz),
                       x = (0, Lx),
                       z = (-Lz, 0),
                       topology=(Bounded, Flat, Bounded))

#################################### Boundary conditions ####################################

# Wind stress boundary condition
# τ₀ = 1e-5           # m² s⁻²
@inline function τy(x, t, p)
    return 0.0125/1e3 * (sinpi(x/p.Lx/0.8-0.5)+1)*(-tanh(x/p.Lx*pi/0.1-9*pi)+1)
end

y_wind_stress = FluxBoundaryCondition(τy, parameters=(; Lx)) #τ₀,
v_bcs = FieldBoundaryConditions(top=y_wind_stress)

#################################### Model ####################################

background_diffusivity = ScalarDiffusivity(ν=1e-4, κ=1e-4)
catke = CATKEVerticalDiffusivity()
# horizontal_closure = HorizontalScalarDiffusivity(ν=1e3, κ=1e3)

# Model
model = HydrostaticFreeSurfaceModel(; grid,
                                    buoyancy = BuoyancyTracer(),
                                    coriolis = FPlane(; f=1e-4),
                                    closure = (catke, background_diffusivity), #horizontal_closure),
                                    tracers = (:b, :e, :T1, :T2),
                                    momentum_advection = WENO(),
                                    tracer_advection = WENO(),
                                    boundary_conditions = (; v=v_bcs))

M² = 0         # 1e-7 s⁻², squared buoyancy frequency
N² = 0         # 1e-5 s⁻²
bᵢ(x, z) = N² * z + M² * x

# Add a layered passive tracer c to track circulation patterns
T1ᵢ(x, z) = floor(Int, -z / 1000) % 2

# hᵢ(x, z) =  floor(Int, x / 500kilometers) % 2

function create_tracer(Nx, Nz) # 100, 200
    # Create the tracer field initialized with zeros
    tracer = zeros(Nx, 1, Nz)

    # Set the tracer value to 1 at the specified locations
    tracer[45:55,1,190:200] .= 1
    tracer[80:90,1,95:105] .= 1

    return tracer
end

set!(model, b=bᵢ, T1 = T1ᵢ, T2 = create_tracer(Nx, Nz)) 

simulation = Simulation(model; Δt = 30minutes, stop_time=3days)

outputs = (u = model.velocities.u,
            w = model.velocities.w,
            T1 = model.tracers.T1,
            T2 = model.tracers.T2)

simulation.output_writers[:simple_output] =
    JLD2OutputWriter(model, outputs, 
                     schedule = TimeInterval(1hour), #30minutes
                     filename = "test",
                     overwrite_existing = true)

run!(simulation)

############################ Visualizing the solution ############################

# filepath = simulation.output_writers[:simple_output].filepath

u_timeseries = FieldTimeSeries(filepath, "u")
w_timeseries = FieldTimeSeries(filepath, "w")
times = w_timeseries.times
xw, yw, zw = nodes(w_timeseries)

T1_timeseries = FieldTimeSeries(filepath, "T1")
xt, yt, zt = nodes(T1_timeseries)
T2_timeseries = FieldTimeSeries(filepath, "T2")

n = Observable(1)
title = @lift @sprintf("t = %d hour", times[$n] / hour) 

uₙ = @lift interior(u_timeseries[$n], :, 1, :)
wₙ = @lift interior(w_timeseries[$n], :, 1, :)
T1ₙ = @lift interior(T1_timeseries[$n], :, 1, :)
T2ₙ = @lift interior(T2_timeseries[$n], :, 1, :)

fig = Figure(size = (1200, 1200))

ax_u = Axis(fig[2, 1]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_u = heatmap!(ax_u, xw, zw, uₙ; colorrange = (-1e-3, 1e-3), colormap = :balance) 
Colorbar(fig[2, 2], hm_u; label = "u", flipaxis = false)

ax_w = Axis(fig[2, 3]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_w = heatmap!(ax_w, xw, zw, wₙ; colorrange = (-2e-6,2e-6), colormap = :balance) 
Colorbar(fig[2, 4], hm_w; label = "w", flipaxis = false)

ax_T1 = Axis(fig[3, 1]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_T1 = heatmap!(ax_T1, xt, zt, T1ₙ; colorrange = (0.4,0.6), colormap = :rainbow1) 
Colorbar(fig[3, 2], hm_T1; label = "Passive tracer (top-down)", flipaxis = false)

ax_T2 = Axis(fig[3, 3]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_T2 = heatmap!(ax_T2, xt, zt, T2ₙ; colorrange = (0, 0.1), colormap = :rainbow1) 
Colorbar(fig[3, 4], hm_T2; label = "Passive tracer (two blocks)", flipaxis = false)

fig[1, 1:4] = Label(fig, title, tellwidth=false)

# And, finally, we record a movie.
frames = 1:length(times)
record(fig, "test.mp4", frames, framerate=10) do i
    n[] = i
end
nothing #hide

# fig2 = Figure(size = (500, 500))
# t2ₙ = @lift interior(T2_timeseries[$n], 38:59, 1, 180:200)
# ax_t2 = Axis(fig2[2, 1]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
# hm_t2 = heatmap!(ax_t2, xt[38:59], zt[180:200], t2ₙ; colorrange = (0,0.02), colormap = :rainbow1) 
# Colorbar(fig2[2, 2], hm_t2; label = "Passive tracer (zoom in)", flipaxis = false)

# fig2[1, 1:2] = Label(fig2, title, tellwidth=false)

# frames = 1:length(times)
# record(fig2, "test.mp4", frames, framerate=50) do i
#     n[] = i
# end
# nothing #hide