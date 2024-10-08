# # test POP remineralization formulations in 2D settings
#
# The BGC tracers in our model are advected, diffuse, precipitate, and dissolve 

using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients
# using GLMakie
using CUDA
using Printf
using Statistics

using Oceananigans
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity, CATKEMixingLength
using Oceananigans.Units

#################################### Grid ####################################
# Parameters
Nx = 100 
Nz = 200
Lx = 1000kilometers   # m
Lz = 2000           # m

# We use a two-dimensional grid, with a `Flat` `y`-direction:
grid = RectilinearGrid(GPU(),
                       size = (Nx, Nz),
                       x = (0, Lx),
                       z = (-Lz, 0),
                       topology=(Bounded, Flat, Bounded))

#################################### Boundary conditions ####################################

const τ₀ = 2e-4           # m² s⁻²
function τy(x, t, p)
    if x/p.Lx < 0.8
        return p.τ₀ * sinpi(x/p.Lx / 1.6)
    elseif x/p.Lx >=0.8
        return p.τ₀ * sinpi((1 - x/p.Lx) / 0.4)
    end
end
y_wind_stress = FluxBoundaryCondition(τy, parameters=(; τ₀, Lx)) 
v_bcs = FieldBoundaryConditions(top=y_wind_stress)

#################################### Model ####################################

# mixing_length = CATKEMixingLength(Cᵇ = 0.001)
# catke = CATKEVerticalDiffusivity(; mixing_length, tke_time_step = 10minutes)

kz(x,z,t) = 5e-4 * (tanh((z+150)/20)+1) + 1e-5
vertical_closure = VerticalScalarDiffusivity(;ν=1e-4, κ=kz)
horizontal_closure = HorizontalScalarDiffusivity(ν=1e3)

# Model
model = HydrostaticFreeSurfaceModel(; grid,
                                    buoyancy = BuoyancyTracer(),
                                    biogeochemistry = CarbonAlkalinityNutrients(; grid,
                                                                        maximum_net_community_production_rate  = 1e-2/day,
                                                                        iron_scavenging_rate = 0),
                                                                        #option_of_particulate_remin = 2),
                                    coriolis = FPlane(; f=1e-4),
                                    closure = (vertical_closure, horizontal_closure), 
                                    tracers = (:b, :e, :DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe),
                                    momentum_advection = WENO(),
                                    tracer_advection = WENO(),
                                    boundary_conditions = (; v=v_bcs)) 

M² = 0         # 1e-7 s⁻², squared buoyancy frequency
N² = 0         # 1e-5 s⁻²
bᵢ(x, z) = N² * z + M² * x

# Biogeochemical
# N₀ = 1e-3 # Surface nutrient concentration = 1 uM
# D₀ = 1e-3 # Surface detritus concentration = 1 uM

# Nᵢ(x,z) = N₀ * (-z / 115).^0.84
# Dᵢ(x,z) = D₀ * (-z / 115).^-0.84

set!(model, b=bᵢ, DIC=2.1, ALK=2.35, NO₃=3e-2, PO₄=2e-3, DOP=0, POP=0, Fe = 1e-6) # mol PO₄ m⁻³

simulation = Simulation(model; Δt = 5minutes, stop_time=365.25*50days) 

# We add a `TimeStepWizard` callback to adapt the simulation's time-step,
wizard = TimeStepWizard(cfl=0.2, max_change=1.1, max_Δt=30minutes)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(50))

# Print the progress 
progress(sim) = @printf("Iteration: %d, time: %s\n", #, total(P): %.2e
                        iteration(sim), prettytime(sim)) 

add_callback!(simulation, progress, IterationInterval(500))

# outputs = merge(model.velocities, model.tracers)
outputs = (u = model.velocities.u,
            v = model.velocities.v,
            w = model.velocities.w,
            PO₄= model.tracers.PO₄,
            DOP = model.tracers.DOP,
            POP = model.tracers.POP,
            NO₃ = model.tracers.NO₃,
            Fe = model.tracers.Fe,
            DIC = model.tracers.DIC,
            ALK = model.tracers.ALK,
            avg_PO₄ = Average(model.tracers.PO₄, dims=(1, 2)),
            avg_DOP = Average(model.tracers.DOP, dims=(1, 2)),
            avg_POP = Average(model.tracers.POP, dims=(1, 2)))

simulation.output_writers[:simple_output] =
    JLD2OutputWriter(model, outputs, 
                     schedule = TimeInterval(7days), 
                     filename = "CAN_2D_Rz2",
                     overwrite_existing = true)

simulation.output_writers[:checkpointer] = Checkpointer(model,
                    schedule = TimeInterval(365.25*10days),
                    prefix = "CAN_2D_Rz_checkpoint",
                    overwrite_existing = true)

run!(simulation, pickup = false)

############################ Visualizing the solution ############################

# filepath = simulation.output_writers[:simple_output].filepath

# u_timeseries = FieldTimeSeries(filepath, "u")
# w_timeseries = FieldTimeSeries(filepath, "w")
# times = w_timeseries.times
# xw, yw, zw = nodes(w_timeseries)

# PO4_timeseries = FieldTimeSeries(filepath, "PO₄")
# avg_PO4_timeseries = FieldTimeSeries(filepath, "avg_PO₄")
# POP_timeseries = FieldTimeSeries(filepath, "POP")
# avg_POP_timeseries = FieldTimeSeries(filepath, "avg_POP")
# xi, yi, zi = nodes(avg_PO4_timeseries)

# DOP_timeseries = FieldTimeSeries(filepath, "DOP")
# DIC_timeseries = FieldTimeSeries(filepath, "DIC")
# ALK_timeseries = FieldTimeSeries(filepath, "ALK")
# Fe_timeseries = FieldTimeSeries(filepath, "Fe")
# NO3_timeseries = FieldTimeSeries(filepath, "NO₃")

# n = Observable(1)
# # title = @lift @sprintf("t = Week %d", times[$n] / 7days) 
# title = @lift @sprintf("t = Day %d", times[$n] / 1day) 

# # convert unit from mol/m³ to μM: 1e3*interior(...)
# PO4ₙ = @lift 1e3*interior(PO4_timeseries[$n], :, 1, :)
# avg_PO4ₙ = @lift 1e3*interior(avg_PO4_timeseries[$n], 1, 1, :)

# POPₙ = @lift 1e3*interior(POP_timeseries[$n], :, 1, :)
# avg_POPₙ = @lift 1e3*interior(avg_POP_timeseries[$n], 1, 1, :)

# # # convert unit from m/s to cm/s:
# uₙ = @lift 100*interior(u_timeseries[$n], :, 1, :)
# wₙ = @lift 100*interior(w_timeseries[$n], :, 1, :)

# fig = Figure(size = (1200, 1500))

# ax_u = Axis(fig[2, 1]; xlabel = "x (m)", ylabel = "z (m)",title = "u (cm/s)", aspect = 1)
# hm_u = heatmap!(ax_u, xw, zw, uₙ; colorrange = (-5e-1, 4e-1), colormap = :balance) 
# Colorbar(fig[2, 2], hm_u; flipaxis = false)

# ax_w = Axis(fig[2, 3]; xlabel = "x (m)", ylabel = "z (m)", title = "w (cm/s)", aspect = 1)
# hm_w = heatmap!(ax_w, xw, zw, wₙ; colorrange = (-2e-3,2e-3), colormap = :balance) 
# Colorbar(fig[2, 4], hm_w; flipaxis = false)

# ax_PO4 = Axis(fig[3, 1]; xlabel = "x (m)", ylabel = "z (m)", title = "[PO₄] (μM)", aspect = 1)
# hm_PO4 = heatmap!(ax_PO4, xw, zw, PO4ₙ; colorrange = (0,4),colormap = :matter) 
# Colorbar(fig[3, 2], hm_PO4; flipaxis = false)

# ax_POP = Axis(fig[3, 3]; xlabel = "x (m)", ylabel = "z (m)", title = "[POP] (μM)", aspect = 1)
# hm_POP = heatmap!(ax_POP, xw, zw, POPₙ; colorrange = (0,0.1),colormap = :matter) 
# Colorbar(fig[3, 4], hm_POP; flipaxis = false)

# ax_avg_PO4 = Axis(fig[4, 1:2]; xlabel = "[PO₄] (μM)", ylabel = "z (m)", title = "Average [PO₄] (μM)", yaxisposition = :right)
# xlims!(ax_avg_PO4, 0, 4)
# lines!(ax_avg_PO4, avg_PO4ₙ, zi)

# ax_avg_POP = Axis(fig[4, 3:4]; xlabel = "[POP] (μM)", ylabel = "z (m)", title = "Average [POP] (μM)",yaxisposition = :right)
# xlims!(ax_avg_POP, 0, 0.05)
# lines!(ax_avg_POP, avg_POPₙ, zi)

# fig[1, 1:4] = Label(fig, title, tellwidth=false)

# # And, finally, we record a movie.
# frames = 1:length(times)
# record(fig, "CAN_2D_Rz2.mp4", frames, framerate=50) do i
#     n[] = i
# end
# nothing #hide

################################## Plot all BGC tracers ##################################

# PO4_timeseries = FieldTimeSeries(filepath, "PO₄")
# POP_timeseries = FieldTimeSeries(filepath, "POP")
# DOP_timeseries = FieldTimeSeries(filepath, "DOP")
# DIC_timeseries = FieldTimeSeries(filepath, "DIC")
# ALK_timeseries = FieldTimeSeries(filepath, "ALK")
# Fe_timeseries = FieldTimeSeries(filepath, "Fe")
# NO3_timeseries = FieldTimeSeries(filepath, "NO₃")

# PO4_final = 1e3*interior(PO4_timeseries[end], :, 1, :)
# POP_final = 1e3*interior(POP_timeseries[end], :, 1, :)
# DOP_final = 1e3*interior(DOP_timeseries[end], :, 1, :)
# DIC_final = 1e3*interior(DIC_timeseries[end], :, 1, :)
# ALK_final = 1e3*interior(ALK_timeseries[end], :, 1, :)
# Fe_final = 1e6*interior(Fe_timeseries[end], :, 1, :)
# NO3_final = 1e3*interior(NO3_timeseries[end], :, 1, :)

#=
fig_can = Figure(size = (1000, 1000))

ax_PO4 = Axis(fig_can[2, 1]; xlabel = "x (m)", ylabel = "z (m)", title = "[PO₄] (μM)", aspect = 1)
hm_PO4 = heatmap!(ax_PO4, xw, zw, PO4_final; colorrange = (0,4),colormap = :matter) 
Colorbar(fig_can[2, 2], hm_PO4; flipaxis = false)

ax_POP = Axis(fig_can[2, 3]; xlabel = "x (m)", ylabel = "z (m)", title = "[POP] (μM)", aspect = 1)
hm_POP = heatmap!(ax_POP, xw, zw, POP_final; colorrange = (0,0.15),colormap = :matter) 
Colorbar(fig_can[2, 4], hm_POP; flipaxis = false)

ax_DOP = Axis(fig_can[2, 5]; xlabel = "x (m)", ylabel = "z (m)", title = "[DOP] (μM)", aspect = 1)
hm_DOP = heatmap!(ax_DOP, xw, zw, DOP_final; colorrange = (0,0.15),colormap = :matter) 
Colorbar(fig_can[2, 6], hm_DOP; flipaxis = false)

ax_NO3 = Axis(fig_can[3, 1]; xlabel = "x (m)", ylabel = "z (m)", title = "[NO₃] (μM)", aspect = 1)
hm_NO3 = heatmap!(ax_NO3, xw, zw, NO3_final; colorrange = (0,100), colormap = :matter) 
Colorbar(fig_can[3, 2], hm_NO3; flipaxis = false)

ax_Fe = Axis(fig_can[3, 3]; xlabel = "x (m)", ylabel = "z (m)", title = "[Fe] (nM)", aspect = 1)
hm_Fe = heatmap!(ax_Fe, xw, zw, Fe_final; colorrange = (0,2.5), colormap = :matter) 
Colorbar(fig_can[3, 4], hm_Fe; flipaxis = false)

ax_DIC = Axis(fig_can[4, 1]; xlabel = "x (m)", ylabel = "z (m)", title = "[DIC] (μM)", aspect = 1)
hm_DIC = heatmap!(ax_DIC, xw, zw, DIC_final; colorrange = (1900,2600), colormap = :matter) 
Colorbar(fig_can[4, 2], hm_DIC; flipaxis = false)

ax_ALK = Axis(fig_can[4, 3]; xlabel = "x (m)", ylabel = "z (m)", title = "[ALK] (μM)", aspect = 1)
hm_ALK = heatmap!(ax_ALK, xw, zw, ALK_final; colorrange = (2300,2400), colormap = :matter) 
Colorbar(fig_can[4, 4], hm_ALK; flipaxis = false)

display(fig_can)
=#
################################## Plot prod/remin rates ##################################

# PO4_timeseries = FieldTimeSeries(filepath, "PO₄")
# PO4_final = interior(PO4_timeseries[end], :, 1, :)
# NO3_timeseries = FieldTimeSeries(filepath, "NO₃")
# NO3_final = interior(NO3_timeseries[end], :, 1, :)
# Fe_timeseries = FieldTimeSeries(filepath, "Fe")
# Fe_final = interior(Fe_timeseries[end], :, 1, :)
# POP_timeseries = FieldTimeSeries(filepath, "POP")
# # POP_final = interior(POP_timeseries[end], :, 1, :)
#=
z_matrix = repeat(zi, 1, 100)
μᵖ= 0.01 # /day
kᴵ=10
kᴾ=1e-7*1024.5
kᴺ=1.6e-6*1024.5
kᶠ=1e-10*1024.5
I =700 .* exp.(z_matrix' ./ 25)
P = PO4_final
N = NO3_final
F = Fe_final

light_lim = I ./ (I .+ kᴵ)
p_lim = (P ./ (P .+ kᴾ))
n_lim = (N ./ (N .+ kᴺ))
f_lim = (F ./ (F .+ kᶠ))

NCP_final = max.(0, μᵖ .* light_lim .* min.(p_lim, n_lim, f_lim))
   
z₀ = log(0.01)*25 
# Remin_final = -0.84 .* 10 ./(z_matrix' .+ z₀) .* POP_final
Remin_final = 0.03 .* POP_final

fig_rate = Figure(size = (1000, 1200))

ax_NCP = Axis(fig_rate[1, 1]; xlabel = "x (m)", ylabel = "z (m)",title = "NCP (mmol m⁻³ d⁻¹)", aspect = 1)
hm_NCP = heatmap!(ax_NCP, xw, zw, 1000*NCP_final; colorrange = (0,10),colormap = :matter) 
Colorbar(fig_rate[1, 2], hm_NCP; flipaxis = false)

ax_remin = Axis(fig_rate[1, 3]; xlabel = "x (m)", ylabel = "z (m)",title = "Remin (mmol m⁻³ d⁻¹)", aspect = 1)
hm_remin = heatmap!(ax_remin, xw, zw, 1000*Remin_final; colorrange = (0,2),colormap = :matter) 
Colorbar(fig_rate[1, 4], hm_remin; flipaxis = false)

ax_logNCP = Axis(fig_rate[2, 1]; xlabel = "x (m)", ylabel = "z (m)",title = "Log(NCP) (mmol m⁻³ d⁻¹)", aspect = 1)
hm_logNCP = heatmap!(ax_logNCP, xw, zw[170:200], log.(abs.(1000*NCP_final[:,170:200])); colorrange = (-2,2.5),colormap = :matter) 
Colorbar(fig_rate[2, 2], hm_logNCP; flipaxis = false)

ax_logremin = Axis(fig_rate[2, 3]; xlabel = "x (m)", ylabel = "z (m)",title = "Log(Remin) (mmol m⁻³ d⁻¹)", aspect = 1)
hm_logremin = heatmap!(ax_logremin, xw, zw, log.(abs.(1000*Remin_final)); colorrange = (-10,1),colormap = :matter) 
Colorbar(fig_rate[2, 4], hm_logremin; flipaxis = false)

avg_NCP = mean(1000*NCP_final; dims = 1)
ax_avg_NCP = Axis(fig_rate[3, 1:2]; xlabel = "Rate (mmol m⁻³ d⁻¹)", ylabel = "z (m)", title = "Ave. NCP", yaxisposition = :right)
xlims!(ax_avg_NCP, -0.001, 6.6)
lines!(ax_avg_NCP, vec(avg_NCP),zi)

avg_R = mean(1000*Remin_final; dims = 1)
ax_avg_R = Axis(fig_rate[3, 3:4]; xlabel = "Rate (mmol m⁻³ d⁻¹)", ylabel = "z (m)", title = "Ave. Remin", yaxisposition = :right)
xlims!(ax_avg_R, 0, 1.1)
lines!(ax_avg_R, vec(avg_R),zi)
display(fig_rate)

# limitation terms
fig_lim = Figure(size = (800, 800))

ax_l_lim = Axis(fig_lim[1, 1]; xlabel = "x (m)", ylabel = "z (m)",title = "I/(I+kᵢ)", aspect = 1)
hm_l_lim = heatmap!(ax_l_lim, xw, zw, light_lim; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_lim[1, 2], hm_l_lim; flipaxis = false)

ax_p_lim = Axis(fig_lim[1, 3]; xlabel = "x (m)", ylabel = "z (m)",title = "PO₄/(PO₄+kₚₒ₄)", aspect = 1)
hm_p_lim = heatmap!(ax_p_lim, xw, zw, p_lim; colorrange = (0, 1.1), colormap = :tempo) 
Colorbar(fig_lim[1, 4], hm_p_lim; flipaxis = false)

ax_n_lim = Axis(fig_lim[2, 1]; xlabel = "x (m)", ylabel = "z (m)",title = "NO₃/(NO₃+kₙₒ₃)", aspect = 1)
hm_n_lim = heatmap!(ax_n_lim, xw, zw, n_lim; colorrange = (0, 1.1), colormap = :tempo) 
Colorbar(fig_lim[2, 2], hm_n_lim; flipaxis = false)

ax_f_lim = Axis(fig_lim[2, 3]; xlabel = "x (m)", ylabel = "z (m)",title = "Fe/(Fe+kᵢᵣₒₙ)", aspect = 1)
hm_f_lim = heatmap!(ax_f_lim, xw, zw, f_lim; colorrange = (0, 1.1), colormap = :tempo) 
Colorbar(fig_lim[2, 4], hm_f_lim; flipaxis = false)

display(fig_lim)

# compare to Martin curve
Remin_flux = vec(mean(10 .* POP_final; dims = 1))
Martin_flux = Remin_flux[175]*((zi[175]+z₀)./(zi.+z₀)).^0.84

fig_flux = Figure(size = (500, 500))
ax_flux = Axis(fig_flux[1, 1]; xlabel = "POP flux (mmol m⁻² d⁻¹)", ylabel = "z (m)", title = "Flux comparison", yaxisposition = :right)
xlims!(ax_flux, 0, 0.5)
lines!(ax_flux, Remin_flux, zi, label = "Model")
lines!(ax_flux, Martin_flux, zi, label = "Martin curve")
axislegend(ax_flux, position = :rb)

display(fig_flux)
=#
################################## Plot the final velocities ##################################

# filepath = "./CAN_NormalWind_NCP1e_2_10y.jld2"

# Plot the final velocity field in log scale: 
# u_final = (interior(u_timeseries[end], :, 1, :))
# w_final = (interior(w_timeseries[end], 1:80, 1, :))

# fig_v = Figure(size = (600, 600))

# ax_uf = Axis(fig_v[1, 1]; xlabel = "x (m)", ylabel = "z (m)",title = "u (m/s)", aspect = 1)
# hm_uf = heatmap!(ax_uf, xw, zw, u_final; colorrange = (-3e-4, 3e-4), colormap = :balance) 
# Colorbar(fig_v[1, 2], hm_uf; flipaxis = false)

# ax_wf = Axis(fig_v[1, 1]; xlabel = "x (m)", ylabel = "z (m)", title = "w (m/s)", aspect = 1)
# hm_wf = heatmap!(ax_wf, xw[1:80], zw, w_final; colorrange = (-1.2e-6, 1.2e-6), colormap = :balance) 
# Colorbar(fig_v[1, 2], hm_wf; flipaxis = false)
# display(fig_v)

################################## Last P cycle ##################################

# PO4_last = 1e3*interior(avg_PO4_timeseries[end], 1, 1, :) 
# POP_last = 1e3*interior(avg_POP_timeseries[end], 1, 1, :) 

# last_frame = Figure(size=(1000, 500))

# ax_avePO4  = Axis(last_frame[1, 1], xlabel="[PO₄] (μM)", ylabel="Depth (m)")
# # xlims!(ax_avePO4, 16,20)
# lines!(ax_avePO4, PO4_last, zi, linewidth = 2)

# ax_avePOP  = Axis(last_frame[1, 2], xlabel="[POP] (μM)", ylabel="Depth (m)")
# xlims!(ax_avePOP, 0, 5e-2)
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
