# AMOC-like simulation

# using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients
using GLMakie
# using CUDA
using Printf
using Statistics

using Oceananigans
# using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity, CATKEMixingLength
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization
using Oceananigans.Units

#################################### Grid ####################################
# Parameters
Nx = 500 
Nz = 100
Lx = 15000kilometers   # m
Lz = 2000           # m

# We use a two-dimensional grid, with a `Flat` `y`-direction:
grid = RectilinearGrid(#GPU(),
                       size = (Nx, Nz),
                       x = (0, Lx),
                       z = (-Lz, 0),
                       topology=(Bounded, Flat, Bounded))

#################################### Boundary conditions ####################################

τ₀ = - 1e-4          # m² s⁻²
function τy(x, t, p)
    xfrac = x / p.Lx
    # return xfrac < 0.8 ?  p.τ₀ *sinpi(5*xfrac/8) : p.τ₀ * cospi((2.5 * xfrac - 2))
    return p.τ₀ * sinpi(xfrac)
end
y_wind_stress = FluxBoundaryCondition(τy, parameters=(; τ₀=τ₀, Lx=Lx)) 
v_bcs = FieldBoundaryConditions(top=y_wind_stress)

#################################### Model ####################################

kz(x,z,t) = 1e-5 + 5e-3 * (tanh((z+150)/20)+1) + 1e-2 * exp(-(z+2000)/50)
vertical_closure = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(),ν=1, κ=kz)
horizontal_closure = HorizontalScalarDiffusivity(ν=1e3, κ=1e3)

# Model
model = HydrostaticFreeSurfaceModel(; grid,
                                    buoyancy = BuoyancyTracer(),
                                    coriolis = FPlane(; f=1e-5),
                                    closure = (vertical_closure, horizontal_closure), 
                                    tracers = (:b),
                                    # momentum_advection = WENO(),#CenteredSecondOrder()
                                    # tracer_advection = WENO(),
                                    boundary_conditions = (; v=v_bcs)) 

M² = 0        # 1e-7 s⁻², squared buoyancy frequency
N² = 0         # 1e-5 s⁻²
bᵢ(x, z) = N² * z + M² * x

set!(model, b=bᵢ) 

simulation = Simulation(model; Δt = 5minutes, stop_time=100days) 

# We add a `TimeStepWizard` callback to adapt the simulation's time-step,
wizard = TimeStepWizard(cfl=0.2, max_change=1.1, max_Δt=30minutes)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(50))

# Print the progress 
progress(sim) = @printf("Iteration: %d, time: %s\n", 
                        iteration(sim), prettytime(sim)) 

add_callback!(simulation, progress, IterationInterval(1000))

# outputs = merge(model.velocities, model.tracers)
outputs = (u = model.velocities.u,
            v = model.velocities.v,
            w = model.velocities.w)

simulation.output_writers[:simple_output] =
    JLD2OutputWriter(model, outputs, 
                     schedule = TimeInterval(1days), 
                     filename = "AMOC_physics",
                     overwrite_existing = true)

# simulation.output_writers[:checkpointer] = Checkpointer(model,
#                     schedule = TimeInterval(365.25*250days),
#                     prefix = "AMOC3_checkpoint",
#                     overwrite_existing = true)

run!(simulation) #, pickup = false)

############################ Visualizing the solution ############################

filepath = simulation.output_writers[:simple_output].filepath
# filepath = "./AMOC6.jld2"

u_timeseries = FieldTimeSeries(filepath, "u")
w_timeseries = FieldTimeSeries(filepath, "w")

times = w_timeseries.times
xw, yw, zw = nodes(w_timeseries)

n = Observable(1)
# title = @lift @sprintf("t = Year %d0", times[$n] / 3652.5days) 
title = @lift @sprintf("t = Day %d", times[$n] / 1days) 

# convert unit from m/s to cm/s:
uₙ = @lift 100*interior(u_timeseries[$n], :, 1, :)
wₙ = @lift 100*interior(w_timeseries[$n], :, 1, :)

top_uₙ = @lift 100*interior(u_timeseries[$n], :, 1, 75:100)
top_wₙ = @lift 100*interior(w_timeseries[$n], :, 1, 75:100)

fig = Figure(size = (1000, 1000))

ax_viscosity = Axis(fig[2, 1]; xlabel = " Vertical ν (m² s⁻¹)", ylabel = "z (m)", aspect = 1)
# function z_viscosity(z)
#     return 1e-5 .+ 5e-3 .* (tanh.((z.+150)./20).+1) .+ 1e-2 .* exp.(-(z.+2000)./50)
# end
lines!(ax_viscosity, ones(size(zw)),zw)
ylims!(ax_viscosity,-2000,0)
# display(fig)

ax_wind = Axis(fig[2, 3]; xlabel = "x (km)", ylabel = "τy (N m⁻²)", aspect = 1)
function y_wind(x,Lx)
    xfrac = x ./ Lx 
    # return xfrac < 0.8 ?  (1.25*xfrac) : ((5 - 5 * xfrac))
    return sinpi.(xfrac) 
end
lines!(ax_wind, xw/1e3, 1000 .* τ₀ .* y_wind.(xw,Lx))

ax_u = Axis(fig[3, 1]; xlabel = "x (km)", ylabel = "z (m)",title = "u (cm/s)", aspect = 1)
hm_u = heatmap!(ax_u, xw./1e3, zw, uₙ; colorrange = (-0.1, 0.1), colormap = :balance) 
Colorbar(fig[3, 2], hm_u; flipaxis = false)

ax_w = Axis(fig[3, 3]; xlabel = "x (km)", ylabel = "z (m)", title = "w (cm/s)", aspect = 1)
hm_w = heatmap!(ax_w, xw./1e3, zw, wₙ; colorrange = (-2e-5,2e-5), colormap = :balance) 
Colorbar(fig[3, 4], hm_w; flipaxis = false)

ax_utop = Axis(fig[4, 1]; xlabel = "x (km)", ylabel = "z (m)",title = "Upper u (cm/s)", aspect = 1)
hm_utop = heatmap!(ax_utop, xw./1e3, zw[75:100], top_uₙ; colorrange = (-1, 1), colormap = :balance) 
Colorbar(fig[4, 2], hm_utop; flipaxis = false)

ax_wtop = Axis(fig[4, 3]; xlabel = "x (km)", ylabel = "z (m)", title = "Upper w (cm/s)", aspect = 1)
hm_wtop = heatmap!(ax_wtop, xw./1e3, zw[75:100], top_wₙ; colorrange = (-2e-5,2e-5), colormap = :balance) 
Colorbar(fig[4, 4], hm_wtop; flipaxis = false)

fig[1, 1:4] = Label(fig, title, tellwidth=false)

# And, finally, we record a movie.
frames = 1:length(times)
record(fig, "AMOC_physics.mp4", frames, framerate=25) do i
    n[] = i
end
nothing #hide


################################## Plot the final velocities ##################################

# v_timeseries = FieldTimeSeries(filepath, "v")
# v_final = (interior(v_timeseries[end], :, 1, :))

# fig_v = Figure(size = (600, 600))

# ax_vf = Axis(fig_v[1, 1]; xlabel = "x (m)", ylabel = "z (m)",title = "v (m/s)", aspect = 1)
# hm_vf = heatmap!(ax_vf, xw, zw, v_final; colormap = :dense) 
# Colorbar(fig_v[1, 2], hm_vf; flipaxis = false)

# display(fig_v)

# filepath = "./CAN_NormalWind_NCP1e_2_10y.jld2"

# Plot the final velocity field in log scale: 
# u_final = 100*(interior(u_timeseries[1026], :, 1, :))
# w_final = 100*(interior(w_timeseries[1026], :, 1, :))

# fig_v = Figure(size = (600, 600))

# ax_uf = Axis(fig_v[1, 1]; xlabel = "x (m)", ylabel = "z (m)",title = "u (cm/s)", aspect = 1)
# hm_uf = heatmap!(ax_uf, xw, zw, u_final; colorrange = (-0.005, 0.005), colormap = :balance) 
# Colorbar(fig_v[1, 2], hm_uf; flipaxis = false)

# ax_wf = Axis(fig_v[2, 1]; xlabel = "x (m)", ylabel = "z (m)", title = "w (cm/s)", aspect = 1)
# hm_wf = heatmap!(ax_wf, xw, zw, w_final; colorrange = (-1e-5, 1e-5), colormap = :balance) 
# Colorbar(fig_v[2, 2], hm_wf; flipaxis = false)
# display(fig_v)
