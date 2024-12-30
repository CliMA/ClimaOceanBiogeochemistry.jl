using GLMakie
using Printf
using Statistics

using Oceananigans
using Oceananigans.Models.HydrostaticFreeSurfaceModels:
                    HydrostaticFreeSurfaceModel,
                    PrescribedVelocityFields
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization
using Oceananigans.Units
using Oceananigans.BoundaryConditions: fill_halo_regions!

Ny = 500 
Nz = 100
Ly = 15000kilometers   # m
Lz = 2000           # m

arch = CPU()
# We use a two-dimensional grid, with a `Flat` `y`-direction:
grid = RectilinearGrid(arch,
                       size = (Ny, Nz),
                       y = (0, Ly),
                       z = (-Lz, 0),
                       topology=(Flat, Bounded, Bounded))

# Set streamfunction
deltaN = 500kilometers   # North downwelling width
deltaS = 1000kilometers  # South upwelling width
deltaZ = 1000     # Vertical asymmetry
Ψᵢ(y, z)  = -8 * ((1 - exp(-y / deltaS)) * (1 - exp(-(Ly - y) / deltaN)) * 
            sinpi(z / Lz) * exp(z/deltaZ))

# Ψᵢ(y, z) = 10*(y/Ly) * sinpi(-z/Lz) * exp((2/Lz)*z) * (1-exp((y-Ly)/(0.2*Ly)))
# Ψᵢ(x, z) = 0.1 * sinpi(x/Lx) * sinpi(-z/Lz) * exp((2/Lz)*z) * exp((2/Lx)*x)
Ψ = Field{Face, Center, Face}(grid)
set!(Ψ, Ψᵢ)
fill_halo_regions!(Ψ, arch)

# Set velocity field from streamfunction
v = YFaceField(grid)
w = ZFaceField(grid)
v .= - ∂z(Ψ)
w .= + ∂y(Ψ)
fill_halo_regions!(v, arch)
fill_halo_regions!(w, arch)

fig = Figure(size = (1000, 1000))

ax_s = Axis(fig[1, 1]; xlabel = "y (km)", ylabel = "z (m)",title = "Stream function (m² s⁻¹)", aspect = 1)
hm_s = heatmap!(ax_s, grid.yᵃᶜᵃ[1:500]./1e3, grid.zᵃᵃᶜ[1:100], Ψ; colormap = :viridis) 
Colorbar(fig[1, 2], hm_s; flipaxis = false)
contour!(ax_s, grid.yᵃᶜᵃ[1:500]./1e3, grid.zᵃᵃᶜ[0:100], Ψ, levels = 10, color = :black)

ax_v = Axis(fig[2, 1]; xlabel = "y (km)", ylabel = "z (m)",title = "v = - ∂z(Ψ) (m s⁻¹)", aspect = 1)
hm_v = heatmap!(ax_v, grid.yᵃᶜᵃ[1:500]./1e3, grid.zᵃᵃᶜ[1:100], v; colorrange = (-2e-2, 2e-2), colormap = :balance) 
Colorbar(fig[2, 2], hm_v; flipaxis = false)

ax_w = Axis(fig[2, 3]; xlabel = "y (km)", ylabel = "z (m)",title = "w = ∂y(Ψ) (m s⁻¹)", aspect = 1)
hm_w = heatmap!(ax_w, grid.yᵃᶜᵃ[1:500]./1e3, grid.zᵃᵃᶜ[1:100], w; colorrange = (-5e-6, 5e-6), colormap = :balance) 
Colorbar(fig[2, 4], hm_w; flipaxis = false)

display(fig)

############################# Model setup ############################# 
#=
kz(x,z,t) = 5e-4 + 5e-3 * (tanh((z+150)/20)+1) + 1e-2 * exp(-(z+2000)/50)
vertical_closure = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), κ=kz)
horizontal_closure = HorizontalScalarDiffusivity(κ=1e3)

# Model
model = HydrostaticFreeSurfaceModel(grid = grid,
                                    tracers = (:T1),
                                    velocities = PrescribedVelocityFields(; u, w),
                                    coriolis = nothing,
                                    buoyancy = nothing,
                                    closure = (vertical_closure, horizontal_closure))

T1ᵢ = zeros(Nx, 1, Nz)
T1ᵢ[400:440, 1, 35:90] .= 1
T1ᵢ[10:40, 1, 35:90] .= 1

set!(model, T1 = T1ᵢ) 

simulation = Simulation(model; Δt = 1days, stop_time=365.25*100days) 

# Print the progress 
# progress(sim) = @printf("Iteration: %d, time: %s\n", 
#                         iteration(sim), prettytime(sim)) 

# add_callback!(simulation, progress, IterationInterval(1000))

outputs = (u = model.velocities.u,
            w = model.velocities.w,
            T1 = model.tracers.T1)

simulation.output_writers[:simple_output] =
    JLD2OutputWriter(model, outputs, 
                     schedule = TimeInterval(365.25days), 
                     filename = "AMOC_physics_test",
                     overwrite_existing = true)

run!(simulation) 

#################################### Visualize ####################################
filepath = simulation.output_writers[:simple_output].filepath
# filepath = "./AMOC6.jld2"

u_timeseries = FieldTimeSeries(filepath, "u")
w_timeseries = FieldTimeSeries(filepath, "w")

times = w_timeseries.times
xw, yw, zw = nodes(w_timeseries)

T1_timeseries = FieldTimeSeries(filepath, "T1")

n = Observable(1)
title = @lift @sprintf("t = Year %d", times[$n] / 365.25days) 

# convert unit from m/s to cm/s:
uₙ = @lift interior(u_timeseries[$n], :, 1, :)
wₙ = @lift interior(w_timeseries[$n], :, 1, :)
T1ₙ = @lift interior(T1_timeseries[$n], :, 1, :)


ax_u = Axis(fig[3, 3]; xlabel = "x (km)", ylabel = "z (m)",title = "u (t) (m s⁻¹)", aspect = 1)
hm_u = heatmap!(ax_u, xw./1e3, zw, uₙ; colorrange = (-0.01, 0.01), colormap = :balance) 
Colorbar(fig[3, 4], hm_u; flipaxis = false)

ax_w = Axis(fig[4, 3]; xlabel = "x (km)", ylabel = "z (m)", title = "w (t) (m s⁻¹)", aspect = 1)
hm_w = heatmap!(ax_w, xw./1e3, zw, wₙ; colorrange = (-1e-6, 1e-6), colormap = :balance) 
Colorbar(fig[4, 4], hm_w; flipaxis = false)

ax_T1 = Axis(fig[2, 3]; xlabel = "x (km)", ylabel = "z (m)", title = "Passive tracer",aspect = 1)
hm_T1 = heatmap!(ax_T1, xw./1e3, zw, T1ₙ; colorrange = (0,0.2), colormap = :rainbow1) 
Colorbar(fig[2, 4], hm_T1; flipaxis = false)

fig[1, 1:4] = Label(fig, title, tellwidth=false)

# And, finally, we record a movie.
frames = 1:length(times)
record(fig, "AMOC_physics_test.mp4", frames, framerate=25) do i
    n[] = i
end
nothing #hide
=#