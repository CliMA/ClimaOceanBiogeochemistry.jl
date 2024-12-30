# # test POP remineralization formulations in 2D settings
#
# The BGC tracers in our model are advected, diffuse, precipitate, and dissolve 

using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients
using GLMakie
# using CUDA
using Printf
using Statistics

using Oceananigans
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization
using Oceananigans.Units

#################################### Grid ####################################
# Parameters
Nx = 500 
Nz = 100
Lx = 10000kilometers   # m
Lz = 2000           # m

# We use a two-dimensional grid, with a `Flat` `y`-direction:
grid = RectilinearGrid(GPU(),
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

# mixing_length = CATKEMixingLength(Cᵇ = 0.001)
# catke = CATKEVerticalDiffusivity(; mixing_length, tke_time_step = 10minutes)

kz(x,z,t) = 1e-5 + 5e-3 * (tanh((z+75)/10)+1) + 1e-2 * exp(-(z+2000)/50)
vertical_closure = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν=1, κ=kz)
horizontal_closure = HorizontalScalarDiffusivity(ν=1e3, κ=1e3)

# Model
model = HydrostaticFreeSurfaceModel(; grid,
                                    buoyancy = BuoyancyTracer(),
                                    biogeochemistry = CarbonAlkalinityNutrients(; grid,
                                                                        maximum_net_community_production_rate  = 1e-5/day,
                                                                        particulate_organic_phosphorus_sedremin_timescale = 0.5 / day,
                                                                        iron_scavenging_rate = 0),
                                    coriolis = FPlane(; f=1e-4),
                                    closure = (vertical_closure, horizontal_closure), 
                                    tracers = (:b, :e, :DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe),
                                    # momentum_advection = WENO(),
                                    # tracer_advection = WENO(),
                                    boundary_conditions = (; v=v_bcs)) 

M² = 0         # 1e-7 s⁻², squared buoyancy frequency
N² = 0         # 1e-5 s⁻²
bᵢ(x, z) = N² * z + M² * x

# Biogeochemical
# N₀ = 1e-3 # Surface nutrient concentration = 1 uM
# D₀ = 1e-3 # Surface detritus concentration = 1 uM

# Nᵢ(x,z) = N₀ * (-z / 115).^0.84
# Dᵢ(x,z) = D₀ * (-z / 115).^-0.84

set!(model, b=bᵢ, DIC=2.1, ALK=2.35, NO₃=2.5e-2, PO₄=1.5e-3, DOP=0, POP=0, Fe = 8e-7) # mol PO₄ m⁻³

simulation = Simulation(model; Δt = 5minutes, stop_time=365.25*100days) 

# We add a `TimeStepWizard` callback to adapt the simulation's time-step,
wizard = TimeStepWizard(cfl=0.2, max_change=1.1, max_Δt=60minutes)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(50))

# Print the progress 
progress(sim) = @printf("Iteration: %d, time: %s\n", #, total(P): %.2e
                        iteration(sim), prettytime(sim)) 

add_callback!(simulation, progress, IterationInterval(1000))

# outputs = merge(model.velocities, model.tracers)
outputs = (u = model.velocities.u,
            v = model.velocities.v,
            w = model.velocities.w,
            PO₄= model.tracers.PO₄,
            DOP = model.tracers.DOP,
            POP = model.tracers.POP,
            NO₃ = model.tracers.NO₃,
            Fe = model.tracers.Fe)
            # DIC = model.tracers.DIC,
            # ALK = model.tracers.ALK)

simulation.output_writers[:simple_output] =
    JLD2OutputWriter(model, outputs, 
                     schedule = TimeInterval(365.25days), 
                     filename = "AMOC6_100y",
                     overwrite_existing = true)

# simulation.output_writers[:checkpointer] = Checkpointer(model,
#                     schedule = TimeInterval(365.25*200days),
#                     prefix = "AMOC6_checkpoint",
#                     overwrite_existing = true)

run!(simulation) #, pickup = false)

############################ Visualizing the solution ############################

# filepath = simulation.output_writers[:simple_output].filepath
filepath = "./AMOC6_100y_v01.jld2"

u_timeseries = FieldTimeSeries(filepath, "u")
w_timeseries = FieldTimeSeries(filepath, "w")

times = w_timeseries.times
xw, yw, zw = nodes(w_timeseries)

PO4_timeseries = FieldTimeSeries(filepath, "PO₄")
POP_timeseries = FieldTimeSeries(filepath, "POP")

DOP_timeseries = FieldTimeSeries(filepath, "DOP")
Fe_timeseries = FieldTimeSeries(filepath, "Fe")
NO3_timeseries = FieldTimeSeries(filepath, "NO₃")

n = Observable(1)
title = @lift @sprintf("t = Year %d", times[$n] / 365.25days) 
# title = @lift @sprintf("t = Day %d0", times[$n] / 10days) 

# convert unit from mol/m³ to μM: 1e3*interior(...)
PO4ₙ = @lift 1e3*interior(PO4_timeseries[$n], :, 1, :)
avg_PO4ₙ = @lift mean(1e3*interior(PO4_timeseries[$n], :, 1, :), dims=1) 

DOPₙ = @lift 1e3*interior(DOP_timeseries[$n], :, 1, :)
avg_DOPₙ = @lift mean(1e3*interior(DOP_timeseries[$n], :, 1, :), dims=1) 

POPₙ = @lift 1e3*interior(POP_timeseries[$n], :, 1, :)
avg_POPₙ = @lift mean(1e3*interior(POP_timeseries[$n], :, 1, :), dims=1) 

# convert unit from m/s to cm/s:
uₙ = @lift 100*interior(u_timeseries[$n], :, 1, :)
wₙ = @lift 100*interior(w_timeseries[$n], :, 1, :)

fig = Figure(size = (1400, 1200))

# ax_wind = Axis(fig[2, 3]; xlabel = "x (km)", ylabel = "τy (N m⁻²)", aspect = 1)
# function y_wind(x,Lx)
#     xfrac = x ./ Lx 
#     return sinpi.(xfrac) 
# end
# lines!(ax_wind, xw/1e3, 1000 .* τ₀ .* y_wind(xw,Lx))

ax_u = Axis(fig[2, 1]; xlabel = "x (km)", ylabel = "z (m)",title = "u (cm/s)", aspect = 1)
hm_u = heatmap!(ax_u, xw./1e3, zw, uₙ; colorrange = (-20, 20), colormap = :balance) 
Colorbar(fig[2, 2], hm_u; flipaxis = false)

ax_w = Axis(fig[2, 3]; xlabel = "x (km)", ylabel = "z (m)", title = "w (cm/s)", aspect = 1)
hm_w = heatmap!(ax_w, xw./1e3, zw, wₙ; colorrange = (-0.01,0.01), colormap = :balance) 
Colorbar(fig[2, 4], hm_w; flipaxis = false)

ax_PO4 = Axis(fig[3, 1]; xlabel = "x (m)", ylabel = "z (m)", title = "[PO₄] (μM)", aspect = 1)
hm_PO4 = heatmap!(ax_PO4, xw, zw, PO4ₙ; colorrange = (1,2),colormap = :jet1) 
Colorbar(fig[3, 2], hm_PO4; flipaxis = false)

ax_POP = Axis(fig[3, 3]; xlabel = "x (m)", ylabel = "z (m)", title = "[POP] (μM)", aspect = 1)
hm_POP = heatmap!(ax_POP, xw, zw, POPₙ; colorrange = (0,0.01),colormap = :jet1) 
Colorbar(fig[3, 4], hm_POP; flipaxis = false)

ax_DOP = Axis(fig[3, 5]; xlabel = "x (m)", ylabel = "z (m)", title = "[DOP] (μM)", aspect = 1)
hm_DOP = heatmap!(ax_DOP, xw, zw, DOPₙ; colorrange = (0,0.08),colormap = :jet1) 
Colorbar(fig[3, 6], hm_DOP; flipaxis = false)

ax_avg_PO4 = Axis(fig[4, 1:2]; xlabel = "[PO₄] (μM)", ylabel = "z (m)", title = "Average [PO₄] (μM)", yaxisposition = :right)
xlims!(ax_avg_PO4, 0, 3)
PO4_prof = lines!(ax_avg_PO4, avg_PO4ₙ[][1, :], zw[1:100])

ax_avg_POP = Axis(fig[4, 3:4]; xlabel = "[POP] (μM)", ylabel = "z (m)", title = "Average [POP] (μM)",yaxisposition = :right)
xlims!(ax_avg_POP, 0, 0.01)
POP_prof = lines!(ax_avg_POP, avg_POPₙ[][1, :], zw[1:100])

ax_avg_DOP = Axis(fig[4, 5:6]; xlabel = "[DOP] (μM)", ylabel = "z (m)", title = "Average [DOP] (μM)",yaxisposition = :right)
xlims!(ax_avg_DOP, 0, 0.06)
DOP_prof = lines!(ax_avg_DOP, avg_DOPₙ[][1, :], zw[1:100])

fig[1, 1:6] = Label(fig, title, tellwidth=false)

# And, finally, we record a movie.
frames = 1:length(times)
record(fig, "AMOC6_100y_v01.mp4", frames, framerate=25) do i
    n[] = i
    PO4_prof[1] = avg_PO4ₙ[][1, :]
    POP_prof[1] = avg_POPₙ[][1, :]
    DOP_prof[1] = avg_DOPₙ[][1, :]
end
nothing #hide

################################## Plot all BGC tracers ##################################
#=
# data = JLD2.load("path/to/yourfile.jld2")
# ave_POP = data["timeseries/avg_POP/1827"]

PO4_final = 1e3*interior(PO4_timeseries[end], :, 1, :)
POP_final = 1e3*interior(POP_timeseries[end], :, 1, :)
DOP_final = 1e3*interior(DOP_timeseries[end], :, 1, :)
Fe_final = 1e6*interior(Fe_timeseries[end], :, 1, :)
NO3_final = 1e3*interior(NO3_timeseries[end], :, 1, :)

z_matrix = repeat(zw[1:100], 1, 500)
μᵖ= 1e-5 # /day
kᴵ=10
kᴾ=1e-7*1024.5 * 1e3 # unit conversion 1e3
kᴺ=1.6e-6*1024.5 * 1e3 # unit conversion 1e3
kᶠ=1e-10*1024.5 * 1e6 # unit conversion 1e6
I = max.(0, 700 .* exp.(z_matrix' ./ 25))
P = max.(0, PO4_final)
N = max.(0, NO3_final)
F = max.(0, Fe_final)

light_lim = I ./ (I .+ kᴵ)
p_lim = (P ./ (P .+ kᴾ))
n_lim = (N ./ (N .+ kᴺ))
f_lim = (F ./ (F .+ kᶠ))

NCP_final = max.(0, μᵖ .* light_lim .* min.(p_lim, n_lim, f_lim))
   
z₀ = log(0.01)*25 
rₛₑ = 0.1 #/day
Remin_final = -0.84 .* 10 ./(z_matrix' .+ z₀) .* POP_final
# Remin_final = ifelse(z < -1990, rₛₑ .* POP_final, 0.1 .* POP_final)

Remin_final[:,1] = rₛₑ .* POP_final[:,1] # sedimentary remin

# tracer concentration and limitation terms
fig_can = Figure(size = (1000, 1000))

# Light limitation
ax_l_lim = Axis(fig_can[1, 1]; xlabel = "x (km)", ylabel = "z (m)",title = "I/(I+kᵢ)", aspect = 1)
hm_l_lim = heatmap!(ax_l_lim, xw/1e3, zw, light_lim; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_can[1, 2], hm_l_lim; flipaxis = false)
# PO4 concentrations
ax_PO4 = Axis(fig_can[1, 3]; xlabel = "x (km)", ylabel = "z (m)", title = "[PO₄] (μM)", aspect = 1)
hm_PO4 = heatmap!(ax_PO4, xw/1e3, zw, PO4_final; colorrange = (0, 2.8),colormap = :jet1) 
Colorbar(fig_can[1, 4], hm_PO4; flipaxis = false)
# PO4 limitation
ax_p_lim = Axis(fig_can[1, 5]; xlabel = "x (km)", ylabel = "z (m)",title = "PO₄/(PO₄+kₚₒ₄)", aspect = 1)
hm_p_lim = heatmap!(ax_p_lim, xw/1e3, zw, p_lim; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_can[1, 6], hm_p_lim; flipaxis = false)
# POP concentrations
ax_POP = Axis(fig_can[2, 1]; xlabel = "x (km)", ylabel = "z (m)", title = "[POP] (μM)", aspect = 1)
hm_POP = heatmap!(ax_POP, xw/1e3, zw, POP_final; colorrange = (0,0.0007),colormap = :jet1) 
Colorbar(fig_can[2, 2], hm_POP; flipaxis = false)
# NO3 concentrations
ax_NO3 = Axis(fig_can[2, 3]; xlabel = "x (km)", ylabel = "z (m)", title = "[NO₃] (μM)", aspect = 1)
hm_NO3 = heatmap!(ax_NO3, xw/1e3, zw, NO3_final; colorrange = (0,40), colormap = :jet1) 
Colorbar(fig_can[2, 4], hm_NO3; flipaxis = false)
# NO3 limitation
ax_n_lim = Axis(fig_can[2, 5]; xlabel = "x (km)", ylabel = "z (m)",title = "NO₃/(NO₃+kₙₒ₃)", aspect = 1)
hm_n_lim = heatmap!(ax_n_lim, xw/1e3, zw, n_lim; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_can[2, 6], hm_n_lim; flipaxis = false)
# DOP concentrations
ax_DOP = Axis(fig_can[3, 1]; xlabel = "x (km)", ylabel = "z (m)", title = "[DOP] (μM)", aspect = 1)
hm_DOP = heatmap!(ax_DOP, xw/1e3, zw, DOP_final; colorrange = (0,0.008),colormap = :jet1) 
Colorbar(fig_can[3, 2], hm_DOP; flipaxis = false)
# Fe concentrations
ax_Fe = Axis(fig_can[3, 3]; xlabel = "x (km)", ylabel = "z (m)", title = "[Fe] (nM)", aspect = 1)
hm_Fe = heatmap!(ax_Fe, xw/1e3, zw, Fe_final; colorrange = (0,1.5), colormap = :jet1) 
Colorbar(fig_can[3, 4], hm_Fe; flipaxis = false)
# Fe limitation
ax_f_lim = Axis(fig_can[3, 5]; xlabel = "x (km)", ylabel = "z (m)",title = "Fe/(Fe+kᵢᵣₒₙ)", aspect = 1)
hm_f_lim = heatmap!(ax_f_lim, xw/1e3, zw, f_lim; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_can[3, 6], hm_f_lim; flipaxis = false)

display(fig_can)

################################## Plot prod/remin rates ##################################

fig_rate = Figure(size = (800, 1000))

ax_NCP = Axis(fig_rate[1, 1]; xlabel = "x (m)", ylabel = "z (m)",title = "NCP (mmol m⁻³ d⁻¹)", aspect = 1)
hm_NCP = heatmap!(ax_NCP, xw, zw, 1000*NCP_final; colorrange = (0,5e-4),colormap = :jet1) 
Colorbar(fig_rate[1, 2], hm_NCP; flipaxis = false)

ax_remin = Axis(fig_rate[1, 3]; xlabel = "x (m)", ylabel = "z (m)",title = "Remin (mmol m⁻³ d⁻¹)", aspect = 1)
hm_remin = heatmap!(ax_remin, xw, zw, Remin_final; colorrange = (0,2e-5),colormap = :jet1) 
Colorbar(fig_rate[1, 4], hm_remin; flipaxis = false)

# ax_logNCP = Axis(fig_rate[2, 1]; xlabel = "x (m)", ylabel = "z (m)",title = "Log(NCP) (mmol m⁻³ d⁻¹)", aspect = 1)
# hm_logNCP = heatmap!(ax_logNCP, xw, zw[170:200], log.(abs.(1000*NCP_final[:,170:200])); colorrange = (-4,1),colormap = :matter) 
# Colorbar(fig_rate[2, 2], hm_logNCP; flipaxis = false)

avg_NCP = mean(1000*NCP_final; dims = 1)
ax_avg_NCP = Axis(fig_rate[2, 1:2]; xlabel = "Rate (mmol m⁻³ d⁻¹)", ylabel = "z (m)", title = "Ave. NCP", yaxisposition = :right)
xlims!(ax_avg_NCP, -1e-5, 4e-4)
lines!(ax_avg_NCP, vec(avg_NCP),zw[1:100])

avg_R = mean(Remin_final; dims = 1)
ax_avg_R = Axis(fig_rate[2, 3:4]; xlabel = "Rate (mmol m⁻³ d⁻¹)", ylabel = "z (m)", title = "Ave. POP Remin", yaxisposition = :right)
xlims!(ax_avg_R, -1e-6, 2e-5)
lines!(ax_avg_R, vec(avg_R),zw[1:100])

zi = zw[1:100]
# compare to Martin curve
Remin_flux = vec(mean(10 .* POP_final; dims = 1))
Martin_flux = Remin_flux[86]*((zi[86]+z₀)./(zi.+z₀)).^0.84

ax_flux = Axis(fig_rate[3, 3:4]; xlabel = "POP flux (mmol m⁻² d⁻¹)", ylabel = "z (m)", title = "Flux comparison", yaxisposition = :right)
xlims!(ax_flux, 0, 0.006)
lines!(ax_flux, Remin_flux[3:100], zi[3:100], label = "Model")
lines!(ax_flux, Martin_flux, zi, label = "Martin curve")
axislegend(ax_flux, position = :rb)

display(fig_rate)
=#
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

################################## Last P cycle ##################################

# PO4_last = 1e3*interior(avg_PO4_timeseries[end], 1, 1, :) 
# POP_last = 1e3*interior(avg_POP_timeseries[end], 1, 1, :) 

# last_frame = Figure(size=(1000, 500))

# ax_avePO4  = Axis(last_frame[1, 1], xlabel="[PO₄] (μM)", ylabel="Depth (m)")
# # xlims!(ax_avePO4, 16,20)
# lines!(ax_avePO4, PO4_last, zi, linewidth = 2)

# ax_avePOP  = Axis(last_frame[1, 2], xlabel="[POP] (μM)", ylabel="Depth (m)")
# # xlims!(ax_avePOP, 0, 5e-2)
# lines!(ax_avePOP, POP_last, zi, linewidth = 2)
# display(last_frame)

#=
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
