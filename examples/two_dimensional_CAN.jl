# # test POP remineralization formulations in 2D settings
#
# The BGC tracers in our model are advected, diffuse, precipitate, and dissolve 


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
Nz = 200
Lx = 1000kilometers   # m
Lz = 2000           # m

# We use a two-dimensional grid, with a `Flat` `y`-direction:
grid = RectilinearGrid(size = (Nx, Nz),
                       x = (0, Lx),
                       z = (-Lz, 0),
                       topology=(Bounded, Flat, Bounded))

#################################### Boundary conditions ####################################

# Compute piston velocity `wp` from restoring time-scale τ★
# Δz = zspacings(grid, Center())
# τ★ = 1days  
# wp = Δz / τ★ # piston velocity = 10 m/day

# Compute buoyancy flux Jᵇ
# @inline function Jᵇ(x, t, b, p)
#     b★ = p.M² * x   # m s⁻², reference buoyancy
#     return p.wp * (b★ - b) # m² s⁻³
# end
# top_buoyancy_bc = FluxBoundaryCondition(Jᵇ, field_dependencies=:b, parameters=(; M², wp))
# b_bcs = FieldBoundaryConditions(top=top_buoyancy_bc)

# Wind stress boundary condition
# τ₀ = 1e-5           # m² s⁻²
# @inline step(x, xc, dx) = (1 + tanh((x - xc) / dx)) / 2 # between 0 and 1
# @inline τy(x, t, p) = p.τ₀ * x/p.Lx #(1-tanh(4*(1-x/p.Lx)))
@inline function τy(x, t, p)
    return 0.0125/1e3 * (sinpi(x/p.Lx/0.8-0.5)+1)*(-tanh(x/p.Lx*pi/0.1-9*pi)+1)
end

# @inline function τy(x, t, p)
#     fracX = x/p.Lx
#     if fracX < 0.9 
#         return p.τ₀ * (0.1 + fracX* 4/9)
#     else
#         return  p.τ₀ * (0.5 + 5*(fracX-0.9))
#     end
# end
y_wind_stress = FluxBoundaryCondition(τy, parameters=(; Lx)) #τ₀,
v_bcs = FieldBoundaryConditions(top=y_wind_stress)

#################################### Model ####################################

background_diffusivity = ScalarDiffusivity(ν=1e-4, κ=1e-4)
catke = CATKEVerticalDiffusivity()
# horizontal_closure = HorizontalScalarDiffusivity(ν=1e3, κ=1e3)

# Model
model = HydrostaticFreeSurfaceModel(; grid,
                                    buoyancy = BuoyancyTracer(),
                                    biogeochemistry = CarbonAlkalinityNutrients(; grid,
                                                                        maximum_net_community_production_rate  = 5e-3/day),
                                                                        # fraction_of_particulate_export = 1),
                                    coriolis = FPlane(; f=1e-4),
                                    closure = (catke, background_diffusivity), #horizontal_closure),
                                    tracers = (:b, :e, :DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe),
                                    momentum_advection = WENO(),
                                    tracer_advection = WENO(),
                                    boundary_conditions = (; v=v_bcs)) #, b=b_bcs))

M² = 0         # 1e-7 s⁻², squared buoyancy frequency
N² = 0         # 1e-5 s⁻²
bᵢ(x, z) = N² * z + M² * x

# Biogeochemical
N₀ = 1e-3 # Surface nutrient concentration = 1 uM
D₀ = 1e-3 # Surface detritus concentration = 1 uM

Nᵢ(x,z) = N₀ * (-z / 115).^0.84
Dᵢ(x,z) = D₀ * (-z / 115).^-0.84

set!(model, b=bᵢ, DIC=2.1, ALK=2.35, NO₃=1e-2, PO₄=Nᵢ, DOP=Dᵢ, POP=Dᵢ, Fe = 1e-5) # mol PO₄ m⁻³

simulation = Simulation(model; Δt = 30minutes, stop_time=365.25*10days)

progress(sim) = @printf("Iteration: %d, time: %s, total(P): %.2e\n",
                        iteration(sim), prettytime(sim), 
                        sum(model.tracers.PO₄)+sum(model.tracers.DOP)+sum(model.tracers.POP))
# progress(sim) = @printf("Iteration: %d, time: %s, u = %.2e, w=%.2e\n",
#                         iteration(sim), prettytime(sim),
#                         maximum(model.velocities.u), maximum(model.velocities.w))

add_callback!(simulation, progress, IterationInterval(200))

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
                     schedule = TimeInterval(1days), #1days
                     filename = "Remin_2D",
                     overwrite_existing = true)

run!(simulation)

############################ Visualizing the solution ############################

filepath = simulation.output_writers[:simple_output].filepath

u_timeseries = FieldTimeSeries(filepath, "u")
w_timeseries = FieldTimeSeries(filepath, "w")
times = w_timeseries.times
xw, yw, zw = nodes(w_timeseries)

PO4_timeseries = FieldTimeSeries(filepath, "PO₄")
avg_PO4_timeseries = FieldTimeSeries(filepath, "avg_PO₄")
POP_timeseries = FieldTimeSeries(filepath, "POP")
avg_POP_timeseries = FieldTimeSeries(filepath, "avg_POP")
xi, yi, zi = nodes(avg_PO4_timeseries)

n = Observable(1)
title = @lift @sprintf("t = %d day", times[$n] / day) 

# convert unit from mol/m³ to μM
PO4ₙ = @lift 1e3*interior(PO4_timeseries[$n], :, 1, :)
avg_PO4ₙ = @lift 1e3*interior(avg_PO4_timeseries[$n], 1, 1, :)

POPₙ = @lift 1e3*interior(POP_timeseries[$n], :, 1, :)
avg_POPₙ = @lift 1e3*interior(avg_POP_timeseries[$n], 1, 1, :)

uₙ = @lift interior(u_timeseries[$n], :, 1, :)
wₙ = @lift interior(w_timeseries[$n], :, 1, :)

fig = Figure(size = (1200, 1500))

ax_u = Axis(fig[2, 1]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_u = heatmap!(ax_u, xw, zw, uₙ; colorrange = (-1e-3, 1e-3), colormap = :balance) 
Colorbar(fig[2, 2], hm_u; label = "u", flipaxis = false)

ax_w = Axis(fig[2, 3]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_w = heatmap!(ax_w, xw, zw, wₙ; colorrange = (-2e-6,2e-6), colormap = :balance) 
Colorbar(fig[2, 4], hm_w; label = "w", flipaxis = false)

ax_PO4 = Axis(fig[3, 1]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_PO4 = heatmap!(ax_PO4, xw, zw, PO4ₙ; colorrange = (0,5),colormap = :hsv)
Colorbar(fig[3, 2], hm_PO4; label = "[PO₄] (μM)", flipaxis = false)

ax_POP = Axis(fig[3, 3]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
hm_POP = heatmap!(ax_POP, xw, zw, POPₙ; colorrange = (1, 2),colormap = :hsv) 
Colorbar(fig[3, 4], hm_POP; label = "[POP] (μM)", flipaxis = false)

ax_avg_PO4 = Axis(fig[4, 1:2]; xlabel = "[PO₄] (μM)", ylabel = "z (m)", yaxisposition = :right)
xlims!(ax_avg_PO4, 0, 5)
lines!(ax_avg_PO4, avg_PO4ₙ, zi)

ax_avg_POP = Axis(fig[4, 3:4]; xlabel = "[POP] (μM)", ylabel = "z (m)", yaxisposition = :right)
xlims!(ax_avg_POP, 1, 2)
lines!(ax_avg_POP, avg_POPₙ, zi)

fig[1, 1:4] = Label(fig, title, tellwidth=false)

# And, finally, we record a movie.
frames = 1:length(times)
record(fig, "Remin2D_rz.mp4", frames, framerate=24) do i
    n[] = i
end
nothing #hide

################################## Plot the last frame ##################################

#=
PO4_last = 1e3*interior(avg_PO4_timeseries[end], 1, 1, :) 
POP_last = 1e3*interior(avg_POP_timeseries[end], 1, 1, :) 

last_frame = Figure(size=(1000, 500))

ax_avePO4  = Axis(last_frame[1, 1], xlabel="[PO₄] (μM)", ylabel="Depth (m)")
xlims!(ax_avePO4, 16,20)
lines!(ax_avePO4, PO4_last, zi, linewidth = 2)

ax_avePOP  = Axis(last_frame[1, 2], xlabel="[POP] (μM)", ylabel="Depth (m)")
xlims!(ax_avePOP, 4.8,5.2)
lines!(ax_avePOP, POP_last, zi, linewidth = 2)
display(last_frame)

# Calculate POP flux
ax_Flux  = Axis(last_frame[1, 1], xlabel="POP flux (mmol m⁻² d⁻¹)", ylabel="Depth (m)")
POP_flux = POP_last * 10 # POP sinking velocity
z₀ = log(0.01)*25 
depth_index = grid.Nz-3
martin = POP_flux[depth_index]*((zi[depth_index].+z₀)./(zi.+z₀)).^0.84
# martin[depth_index+1:grid.Nz].=NaN

lines!(ax_Flux, POP_flux, zi, linewidth = 2, label = "Modeled POP flux")
lines!(ax_Flux, martin, zi, linewidth = 2,linestyle=:dash, label = "Martin curve")
axislegend(ax_Flux, position = :rb)
# xlims!(ax_Flux, 5.58, 5.62)

# Calculate fractional remin
ax_Frac  = Axis(last_frame[1, 2], xlabel="POP flux at z vs at z₀ (mmol m⁻² d⁻¹)", ylabel="Depth (m)")
FPOP = POP_flux /POP_flux[depth_index] # POP sinking velocity

Fmartin = ((zi[depth_index].+z₀)./(zi.+z₀)).^0.84
# Fmartin[depth_index+1:grid.Nz].=NaN

lines!(ax_Frac, FPOP, zi, linewidth = 2, label = "Fₚₒₚ(z)/Fₚₒₚ(z₀)")
lines!(ax_Frac, Fmartin, zi, linewidth = 2,linestyle=:dash, label = "Martin curve")
axislegend(ax_Frac, position = :rb)

display(last_frame)

# save("po4_dop_pop.png", last_frame)
# nothing #hide

=#