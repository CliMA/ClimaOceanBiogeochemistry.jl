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

# Add Lagrangian Particles
n_particles = 1
x₀ = Lx * rand(n_particles)
y₀ = ones(n_particles)
z₀ = -Lz/2 * ones(n_particles)

using StructArrays
struct Lagrangian{T, P}
    x :: T
    y :: T
    z :: T
    Particles_2D :: P
end
Particles_init = ones(n_particles)
particles = StructArray{Lagrangian}((x₀, y₀, z₀, Particles_init))
lagrangian_particles = LagrangianParticles(particles; tracked_fields=(; Particles_2D))

# Model
model = HydrostaticFreeSurfaceModel(; grid,
                                    buoyancy = BuoyancyTracer(),
                                    coriolis = FPlane(; f=1e-4),
                                    closure = (catke, background_diffusivity), #horizontal_closure),
                                    tracers = (:b, :e, :c, :h),
                                    momentum_advection = WENO(),
                                    particles = lagrangian_particles,
                                    boundary_conditions = (; v=v_bcs))

M² = 0         # 1e-7 s⁻², squared buoyancy frequency
N² = 0         # 1e-5 s⁻²
bᵢ(x, z) = N² * z + M² * x


# Add a layered passive tracer c to track circulation patterns
cᵢ(x, z) = floor(Int, -z / 1000) % 2

# hᵢ(x, z) =  floor(Int, x / 500kilometers) % 2

function create_tracer(Nx, Nz) # 100, 200
    # Create the tracer field initialized with zeros
    tracer = zeros(Nx,1, Nz)

    # Set the tracer value to 1 at the specified locations
    tracer[45:55,1,95:105] .= 1
    tracer[80:90,1,95:105] .= 1

    return tracer
end

set!(model, b=bᵢ, c = cᵢ, h = create_tracer(Nx, Nz))

simulation = Simulation(model; Δt = 30minutes, stop_time=1days)

outputs = (u = model.velocities.u,
            w = model.velocities.w,
            c = model.tracers.c,
            h = model.tracers.h,
            p = model.particles)

simulation.output_writers[:simple_output] =
    JLD2OutputWriter(model, outputs, 
                     schedule = TimeInterval(30minutes), #1days
                     filename = "PassiveTracer_2D",
                     overwrite_existing = true)

run!(simulation)

############################ Visualizing the solution ############################

filepath = simulation.output_writers[:simple_output].filepath

u_timeseries = FieldTimeSeries(filepath, "u")
w_timeseries = FieldTimeSeries(filepath, "w")
times = w_timeseries.times
xw, yw, zw = nodes(w_timeseries)

c_timeseries = FieldTimeSeries(filepath, "c")
xc, yc, zc = nodes(c_timeseries)
h_timeseries = FieldTimeSeries(filepath, "h")

px_timeseries = FieldTimeSeries(filepath, "px")

n = Observable(1)
title = @lift @sprintf("t = %d day", times[$n] / day) 

uₙ = @lift interior(u_timeseries[$n], :, 1, :)
wₙ = @lift interior(w_timeseries[$n], :, 1, :)
cₙ = @lift interior(c_timeseries[$n], :, 1, :)
hₙ = @lift interior(h_timeseries[$n], :, 1, :)

fig = Figure(size = (1500, 1200))

ax_u = Axis(fig[2, 1]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_u = heatmap!(ax_u, xw, zw, uₙ; colorrange = (-1e-3, 1e-3), colormap = :balance) 
Colorbar(fig[2, 2], hm_u; label = "u", flipaxis = false)

ax_w = Axis(fig[2, 3]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_w = heatmap!(ax_w, xw, zw, wₙ; colorrange = (-2e-6,2e-6), colormap = :balance) 
Colorbar(fig[2, 4], hm_w; label = "w", flipaxis = false)

ax_c = Axis(fig[3, 3]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_c = heatmap!(ax_c, xc, zc, cₙ; colorrange = (0.4,0.6), colormap = :rainbow1) 
Colorbar(fig[3, 4], hm_c; label = "Passive tracer (top-down)", flipaxis = false)

ax_h = Axis(fig[3, 1]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_h = heatmap!(ax_h, xc, zc, hₙ; colorrange = (0, 0.1), colormap = :rainbow1) 
Colorbar(fig[3, 2], hm_h; label = "Passive tracer (two blocks)", flipaxis = false)

ax_p = Axis(fig[3, 5]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_p = heatmap!(ax_p, xp, zp, pₙ) 
# Colorbar(fig[3, 6], hm_p; label = "Lagrangian particles", flipaxis = false)

fig[1, 1:5] = Label(fig, title, tellwidth=false)

# And, finally, we record a movie.
frames = 1:length(times)
record(fig, "PassiveTracer_2D.mp4", frames, framerate=80) do i
    n[] = i
end
nothing #hide

