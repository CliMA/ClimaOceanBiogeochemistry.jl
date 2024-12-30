# # test a 2D physical setting with passive tracers

using GLMakie
# using CUDA
using Printf
using Statistics
# using StaticArrays

using Oceananigans
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity, CATKEMixingLength
using Oceananigans.Units

#################################### Grid ####################################
# Parameters
Nx = 500 
Nz = 100
Lx = 10000kilometers   # m
Lz = 2000           # m

# We use a two-dimensional grid, with a `Flat` `y`-direction:
grid = RectilinearGrid(#GPU(),
                       size = (Nx, Nz),
                       x = (0, Lx),
                       z = (-Lz, 0),
                       topology=(Bounded, Flat, Bounded))

#################################### Boundary conditions ####################################

# Wind stress boundary condition
τ₀ = - 1e-4          # m² s⁻²
function τy(x, t, p)
    xfrac = x / p.Lx
    # return xfrac < 0.8 ?  p.τ₀ *sinpi(5*xfrac/8) : p.τ₀ * cospi((2.5 * xfrac - 2))
    return p.τ₀ * sinpi(xfrac)
end
y_wind_stress = FluxBoundaryCondition(τy, parameters=(; τ₀=τ₀, Lx=Lx)) 
v_bcs = FieldBoundaryConditions(top=y_wind_stress)

#################################### Model ####################################

# mixing_length = CATKEMixingLength(Cᵇ = 0.001)
# catke = CATKEVerticalDiffusivity(; mixing_length, tke_time_step = 10minutes)

kz(x,z,t) = 1e-5 + 5e-3 * (tanh((z+150)/20)+1) + 1e-2 * exp(-(z+2000)/50)
vertical_closure = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(),ν=1, κ=kz)
horizontal_closure = HorizontalScalarDiffusivity(ν=1e3, κ=1e3)

# Add Lagrangian Particles
# n_particles = 10
# x₀ = Lx * rand(n_particles) #Lx/2 * ones(n_particles) #* collect([0.1,0.5,0.9]) #rand(n_particles)
# y₀ = ones(n_particles)
# z₀ = -Lz * rand(n_particles) 
# lagrangian_particles = LagrangianParticles(x=x₀, y=y₀, z=z₀) 

# Model
model = HydrostaticFreeSurfaceModel(; grid,
                                    buoyancy = BuoyancyTracer(),
                                    coriolis = FPlane(; f=1e-4),
                                    closure = (vertical_closure, horizontal_closure),
                                    tracers = (:b, :T1, :T2),
                                    # momentum_advection = WENO(),
                                    # tracer_advection = WENO(),
                                    # particles = lagrangian_particles,
                                    boundary_conditions = (; v=v_bcs))

M² = 0         # 1e-7 s⁻², squared buoyancy frequency
N² = 0         # 1e-5 s⁻²
bᵢ(x, z) = N² * z + M² * x

# Add a layered passive tracer c to track circulation patterns
# T1ᵢ(x, z) = floor(Int, -z / 400) % 2
# T2ᵢ(x, z) =  floor(Int, x / 200kilometers) % 2

T1ᵢ = zeros(Nx, 1, Nz)
T1ᵢ[460:470, 1, 80:90] .= 1
T1ᵢ[20:30, 1, 80:90] .= 1

T2ᵢ = zeros(Nx, 1, Nz)
T2ᵢ[240:260, 1, 70:90] .= 1

set!(model, b=bᵢ, T1 = T1ᵢ, T2 = T2ᵢ) 

simulation = Simulation(model; Δt = 5minutes, stop_time=3650days)

wizard = TimeStepWizard(cfl=0.2, max_change=1.1, max_Δt=30minutes)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(50))

progress(sim) = @printf("Iteration: %d, time: %s\n", 
                        iteration(sim), prettytime(sim)) 
add_callback!(simulation, progress, IterationInterval(1000))

outputs = (u = model.velocities.u,
            w = model.velocities.w,
            T1 = model.tracers.T1,
            T2 = model.tracers.T2)
            # p = model.particles)

simulation.output_writers[:simple_output] =
    JLD2OutputWriter(model, outputs, 
                     schedule = TimeInterval(10days), 
                     filename = "passiveT_AMOCν1",
                     overwrite_existing = true)

run!(simulation)

############################ Visualizing the solution ############################

filepath = simulation.output_writers[:simple_output].filepath
# filepath = "./passiveT_ty02_20y_btm1e2diff.jld2"

# TODO: How to record a video of particle movement
# using JLD2
# myfile = jldopen("passiveTP_50y.jld2")
# xₙ = myfile["timeseries/p/1"]

u_timeseries = FieldTimeSeries(filepath, "u")
w_timeseries = FieldTimeSeries(filepath, "w")
times = w_timeseries.times
xw, yw, zw = nodes(w_timeseries)

T1_timeseries = FieldTimeSeries(filepath, "T1")
xt, yt, zt = nodes(T1_timeseries)
T2_timeseries = FieldTimeSeries(filepath, "T2")

n = Observable(1)
title = @lift @sprintf("t = Day %d0", times[$n] / 10days) 

uₙ = @lift 100*interior(u_timeseries[$n], :, 1, :)
wₙ = @lift 100*interior(w_timeseries[$n], :, 1, :)
T1ₙ = @lift interior(T1_timeseries[$n], :, 1, :)
T2ₙ = @lift interior(T2_timeseries[$n], :, 1, :)

fig = Figure(size = (1000, 1000))

# ax_wind = Axis(fig[2, 3]; xlabel = "x (m)", ylabel = "τy (N m⁻²)", aspect = 1)
# function y_wind(x,Lx)
#     xfrac = x / Lx
#     return xfrac < 0.8 ?  τ₀ * sinpi(xfrac / 1.6) :  τ₀ * sinpi((1 - xfrac) / 0.4)
# end
# lines!(ax_wind, xw, 1000 .* τ₀ .* y_wind.(xw,Lx))

ax_u = Axis(fig[2, 1]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_u = heatmap!(ax_u, xw, zw, uₙ; colorrange = (-0.5, 0.5), colormap = :balance) 
Colorbar(fig[2, 2], hm_u; label = "u (cm/s)", flipaxis = false)

ax_w = Axis(fig[2, 3]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_w = heatmap!(ax_w, xw, zw, wₙ; colorrange = (-2e-5,2e-5), colormap = :balance) 
Colorbar(fig[2, 4], hm_w; label = "w (cm/s)", flipaxis = false)

ax_T1 = Axis(fig[3, 1]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_T1 = heatmap!(ax_T1, xt, zt, T1ₙ; colorrange = (0,0.1), colormap = :rainbow1) 
Colorbar(fig[3, 2], hm_T1; label = "Passive tracer 1", flipaxis = false)

ax_T2 = Axis(fig[3, 3]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_T2 = heatmap!(ax_T2, xt, zt, T2ₙ; colorrange = (0,0.1), colormap = :rainbow1) 
Colorbar(fig[3, 4], hm_T2; label = "Passive tracer 2", flipaxis = false)

fig[1, 1:4] = Label(fig, title, tellwidth=false)

# And, finally, we record a movie.
frames = 1:length(times)
record(fig, "passiveT_AMOCν1.mp4", frames, framerate=50) do i
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
# record(fig2, "test2.mp4", frames, framerate=50) do i
#     n[] = i
# end
# nothing #hide
