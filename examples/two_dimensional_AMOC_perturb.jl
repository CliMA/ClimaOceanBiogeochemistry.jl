# # test POP remineralization formulations in 2D settings
#
# The BGC tracers in our model are advected, diffuse, precipitate, and dissolve 

using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients
using GLMakie
# using CUDA
using Printf
using Statistics

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

τ₀ = -1e-4           # m² s⁻²
function τy(x, t, p)
    xfrac = x / p.Lx
    #return xfrac < 0.8 ?  p.τ₀ * sinpi(xfrac / 1.6) :  p.τ₀ * sinpi((1 - xfrac) / 0.4)
    return p.τ₀ * sinpi(xfrac)
end
y_wind_stress = FluxBoundaryCondition(τy, parameters=(; τ₀=τ₀, Lx=Lx)) 
v_bcs = FieldBoundaryConditions(top=y_wind_stress)

#################################### Model ####################################

# mixing_length = CATKEMixingLength(Cᵇ = 0.001)
# catke = CATKEVerticalDiffusivity(; mixing_length, tke_time_step = 10minutes)

kz(x,z,t) = 1e-5 + 5e-3 * (tanh((z+75)/10)+1) + 1e-2 * exp(-(z+2000)/50)
vertical_closure = VerticalScalarDiffusivity(;ν=kz, κ=kz)
horizontal_closure = HorizontalScalarDiffusivity(ν=1e3, κ=1e3)

############################ Reset nutrient concentrations ############################ 

# using JLD2
# data = JLD2.load("AMOC4.jld2")
# pop_keys = filter(key -> startswith(key, "timeseries/POP/") && occursin(r"\d+$", key), keys(data))
# numbers = [parse(Int, match(r"\d+$", key).match) for key in pop_keys]
# max_number = maximum(numbers)

# mid_POP = data["timeseries/POP/$max_number"]
# mid_DOP = data["timeseries/DOP/$max_number"]
# mid_PO₄ = data["timeseries/PO₄/$max_number"]
# mid_NO₃ = data["timeseries/NO₃/$max_number"]
# mid_Fe = data["timeseries/Fe/$max_number"]
# mid_DIC = data["timeseries/DIC/$max_number"]
# mid_ALK = data["timeseries/ALK/$max_number"]

# mid_u = data["timeseries/u/$max_number"]
# mid_v = data["timeseries/v/$max_number"]
# mid_w = data["timeseries/w/$max_number"]

# Set deep nutrient concentrations 
# mid_PO₄[95:100, 1, 195:200] .= 5e-3

M² = 0         # 1e-7 s⁻², squared buoyancy frequency
N² = 0         # 1e-5 s⁻²
bᵢ(x, z) = N² * z + M² * x

model2 = HydrostaticFreeSurfaceModel(; grid,
                                    buoyancy = BuoyancyTracer(),
                                    biogeochemistry = CarbonAlkalinityNutrients(; grid,
                                                                        maximum_net_community_production_rate  = 1e-5/day,
                                                                        particulate_organic_phosphorus_sedremin_timescale = 1 / day,
                                                                        incident_PAR  = 1400,
                                                                        iron_scavenging_rate = 0),
                                    coriolis = FPlane(; f=1e-4),
                                    closure = (vertical_closure, horizontal_closure), 
                                    tracers = (:b, :e, :DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe),
                                    momentum_advection = WENO(),
                                    tracer_advection = WENO(),
                                    boundary_conditions = (; v=v_bcs)) 

set!(model2, b=bᵢ, DIC=2.1, ALK=2.35, NO₃=2e-2, PO₄=1.5e-3, DOP=0, POP=0, Fe = 8e-7) 

simulation2 = Simulation(model2; Δt = 5minutes, stop_time=822000days) 
# We add a `TimeStepWizard` callback to adapt the simulation's time-step,
wizard = TimeStepWizard(cfl=0.2, max_change=1.1, max_Δt=60minutes)
simulation2.callbacks[:wizard] = Callback(wizard, IterationInterval(50))

progress(sim) = @printf("Iteration: %d, time: %s, total(P)=%.2e\n", 
                        iteration(sim), prettytime(sim),                    
                        sum(model2.tracers.PO₄)+sum(model2.tracers.DOP)+sum(model2.tracers.POP)) 
add_callback!(simulation2, progress, IterationInterval(1000))

outputs2 = (u = model2.velocities.u,
            v = model2.velocities.v,
            w = model2.velocities.w,
            PO₄= model2.tracers.PO₄,
            DOP = model2.tracers.DOP,
            POP = model2.tracers.POP,
            NO₃ = model2.tracers.NO₃,
            Fe = model2.tracers.Fe)

simulation2.output_writers[:simple_output] =
    JLD2OutputWriter(model2, outputs2, 
                     schedule = TimeInterval(1days), 
                     filename = "AMOC4_perturbPAR",
                     overwrite_existing = true)

simulation2.output_writers[:checkpointer] = Checkpointer(model2,
                                schedule = TimeInterval(180days),
                                prefix = "AMOC4_checkpoint",
                                overwrite_existing = false)

run!(simulation2, pickup = true)
#################################################################################### 
############################ Visualizing the solution ############################
#################################################################################### 

# filepath = simulation2.output_writers[:simple_output].filepath
filepath = "./AMOC4_perturbPAR.jld2"

u_timeseries = FieldTimeSeries(filepath, "u")
w_timeseries = FieldTimeSeries(filepath, "w")

PO4_timeseries = FieldTimeSeries(filepath, "PO₄")
POP_timeseries = FieldTimeSeries(filepath, "POP")

times = PO4_timeseries.times
xw, yw, zw = nodes(PO4_timeseries)

DOP_timeseries = FieldTimeSeries(filepath, "DOP")
Fe_timeseries = FieldTimeSeries(filepath, "Fe")
NO3_timeseries = FieldTimeSeries(filepath, "NO₃")

n = Observable(1)
title = @lift @sprintf("t = Day %d", times[$n] / 1day) 

# convert unit from mol/m³ to μM: 1e3*interior(...)
PO4ₙ = @lift 1e3*interior(PO4_timeseries[$n], :, 1, :)
DOPₙ = @lift 1e3*interior(DOP_timeseries[$n], :, 1, :)
POPₙ = @lift 1e3*interior(POP_timeseries[$n], :, 1, :)
avg_PO4ₙ = @lift mean(1e3*interior(PO4_timeseries[$n], :, 1, :), dims=1) 
avg_POPₙ = @lift mean(1e3*interior(POP_timeseries[$n], :, 1, :), dims=1) 

# # convert unit from m/s to cm/s:
uₙ = @lift 100*interior(u_timeseries[$n], :, 1, :)
wₙ = @lift 100*interior(w_timeseries[$n], :, 1, :)

fig = Figure(size = (1200, 1500))

ax_u = Axis(fig[2, 1]; xlabel = "x (m)", ylabel = "z (m)",title = "u (cm/s)", aspect = 1)
hm_u = heatmap!(ax_u, xw, zw, uₙ; colorrange = (-0.03, 0.03), colormap = :balance) 
Colorbar(fig[2, 2], hm_u; flipaxis = false)

ax_w = Axis(fig[2, 3]; xlabel = "x (m)", ylabel = "z (m)", title = "w (cm/s)", aspect = 1)
hm_w = heatmap!(ax_w, xw, zw, wₙ; colorrange = (-1e-5,1e-5), colormap = :balance) 
Colorbar(fig[2, 4], hm_w; flipaxis = false)

ax_PO4 = Axis(fig[3, 1]; xlabel = "x (m)", ylabel = "z (m)", title = "[PO₄] (μM)", aspect = 1)
hm_PO4 = heatmap!(ax_PO4, xw, zw, PO4ₙ; colorrange = (0,2.8), colormap = :jet1) 
Colorbar(fig[3, 2], hm_PO4; flipaxis = false)

ax_POP = Axis(fig[3, 3]; xlabel = "x (m)", ylabel = "z (m)", title = "[POP] (μM)", aspect = 1)
hm_POP = heatmap!(ax_POP, xw, zw, POPₙ; colorrange = (0,0.001),colormap = :jet1) 
Colorbar(fig[3, 4], hm_POP; flipaxis = false)

ax_avg_PO4 = Axis(fig[4, 1:2]; xlabel = "[PO₄] (μM)", ylabel = "z (m)", title = "Average [PO₄] (μM)", yaxisposition = :right)
xlims!(ax_avg_PO4, 0, 3)
PO4_prof = lines!(ax_avg_PO4, avg_PO4ₙ[][1, :], zw[1:100])

ax_avg_POP = Axis(fig[4, 3:4]; xlabel = "[POP] (μM)", ylabel = "z (m)", title = "Average [POP] (μM)",yaxisposition = :right)
xlims!(ax_avg_POP, 0, 0.0007)
POP_prof = lines!(ax_avg_POP, avg_POPₙ[][1, :], zw[1:100])

fig[1, 1:4] = Label(fig, title, tellwidth=false)

# And, finally, we record a movie.
frames = 1:length(times)
record(fig, "AMOC4_perturbPAR.mp4", frames, framerate=50) do i
    n[] = i
    PO4_prof[1] = avg_PO4ₙ[][1, :]
    POP_prof[1] = avg_POPₙ[][1, :]
end
nothing #hide

############################## Compare initial vs final ##############################

# PO4_last = 1e3*interior(PO4_timeseries[end], :, 1, :) 
# POP_last = 1e3*interior(POP_timeseries[end], :, 1, :) 
# DOP_last = 1e3*interior(DOP_timeseries[end], :, 1, :) 
# sum(PO4_last.+POP_last.+DOP_last)

PO4_init = 1e3*interior(PO4_timeseries[1], :, 1, :) 
POP_init = 1e3*interior(POP_timeseries[1], :, 1, :) 
# DOP_init = 1e3*interior(DOP_timeseries[1], :, 1, :)
# sum(PO4_init.+POP_init.+DOP_init)

n = Observable(1)
title = @lift @sprintf("t = Day %d", times[$n] / 1days) 

diff_PO4 = @lift 1e3*(interior(PO4_timeseries[$n], :, 1, :) .- interior(PO4_timeseries[1], :, 1, :))
diff_POP = @lift 1e3*(interior(POP_timeseries[$n], :, 1, :) .- interior(POP_timeseries[1], :, 1, :))
# diff_DOP = @lift 1e3*(interior(DOP_timeseries[$n], :, 1, :) .- interior(DOP_timeseries[1], :, 1, :))
# sum(diff_PO4)

# One upwelling station
up_PO4ₙ = @lift 1e3*interior(PO4_timeseries[$n], 100, 1, :)
# up_DOPₙ = @lift 1e3*interior(DOP_timeseries[$n], 100, 1, :)
up_POPₙ = @lift 1e3*interior(POP_timeseries[$n], 100, 1, :)

init_up_PO4ₙ = 1e3*interior(PO4_timeseries[1], 100, 1, :)
init_up_POPₙ = 1e3*interior(POP_timeseries[1], 100, 1, :)
# init_up_DOPₙ = 1e3*interior(DOP_timeseries[1], 100, 1, :)

# One downwelling station
down_PO4ₙ = @lift 1e3*interior(PO4_timeseries[$n], 440, 1, :)
# down_DOPₙ = @lift 1e3*interior(DOP_timeseries[$n], 440, 1, :)
down_POPₙ = @lift 1e3*interior(POP_timeseries[$n], 440, 1, :)

init_down_PO4ₙ = 1e3*interior(PO4_timeseries[1], 440, 1, :)
init_down_POPₙ = 1e3*interior(POP_timeseries[1], 440, 1, :)
# init_down_DOPₙ = 1e3*interior(DOP_timeseries[1], 440, 1, :)

fig_compare = Figure(size=(900, 1100))

# Row 1: final - initial
ax_PO4 = Axis(fig_compare[2, 1]; xlabel = "x (m)", ylabel = "z (m)", title = "Δ[PO₄] (μM)", aspect = 1)
hm_PO4 = heatmap!(ax_PO4, xw, zw, diff_PO4; colorrange = (-0.02,0.02),colormap = :balance) 
Colorbar(fig_compare[2, 2], hm_PO4; flipaxis = false)

ax_POP = Axis(fig_compare[2, 3]; xlabel = "x (m)", ylabel = "z (m)", title = "Δ[POP] (μM)", aspect = 1)
hm_POP = heatmap!(ax_POP, xw, zw, diff_POP; colorrange = (-0.0005,0.0005),colormap = :balance) 
Colorbar(fig_compare[2, 4], hm_POP; flipaxis = false)

# ax_DOP = Axis(fig_compare[2, 5]; xlabel = "x (m)", ylabel = "z (m)", title = "Δ[DOP] (μM)", aspect = 1)
# hm_DOP = heatmap!(ax_DOP, xw, zw, diff_DOP; colorrange = (-0.06,0.06),colormap = :balance) 
# Colorbar(fig_compare[2, 6], hm_DOP; flipaxis = false)

# Row 2: upwelling site
up_PO4 = Axis(fig_compare[3, 1:2]; xlabel = "[PO₄] (μM)", ylabel = "z (m)", yaxisposition = :left)
xlims!(up_PO4, 0, 3)
lines!(up_PO4, init_up_PO4ₙ, zw[1:100],color = :red, linestyle = :dash, label = "Initial")
PO4_up = lines!(up_PO4, up_PO4ₙ[][:, 1], zw[1:100], label = "Dynamic")
axislegend(up_PO4, position = :lb)

up_POP = Axis(fig_compare[3, 3:4]; xlabel = "[POP] (μM)", ylabel = "z (m)", yaxisposition = :left)
xlims!(up_POP, 0, 0.001)
lines!(up_POP, init_up_POPₙ, zw[1:100],color = :red, linestyle = :dash, label = "Initial")
POP_up = lines!(up_POP, up_POPₙ[][:, 1], zw[1:100], label = "Dynamic")
axislegend(up_POP, position = :rb)

# up_DOP = Axis(fig_compare[3, 5:6]; xlabel = "[DOP] (μM)", ylabel = "z (m)", yaxisposition = :left)
# xlims!(up_DOP, -0.02, 0.5)
# lines!(up_DOP, init_up_DOPₙ, zw[1:200],color = :red, linestyle = :dash, label = "Initial")
# DOP_up = lines!(up_DOP, up_DOPₙ[][:, 1], zw[1:200], label = "Dynamic")
# axislegend(up_DOP, position = :rb)

# Row 3: downwelling site
down_PO4 = Axis(fig_compare[4, 1:2]; xlabel = "[PO₄] (μM)", ylabel = "z (m)", yaxisposition = :left)
xlims!(down_PO4, 0, 3)
lines!(down_PO4, init_down_PO4ₙ, zw[1:100],color = :red, linestyle = :dash, label = "Initial")
PO4_down = lines!(down_PO4, down_PO4ₙ[][:, 1], zw[1:100], label = "Dynamic")
axislegend(down_PO4, position = :lb)

down_POP = Axis(fig_compare[4, 3:4]; xlabel = "[POP] (μM)", ylabel = "z (m)", yaxisposition = :left)
xlims!(down_POP, 0, 0.001)
lines!(down_POP, init_down_POPₙ, zw[1:100],color = :red, linestyle = :dash, label = "Initial")
POP_down = lines!(down_POP, down_POPₙ[][:, 1], zw[1:100], label = "Dynamic")
axislegend(down_POP, position = :rb)

# down_DOP = Axis(fig_compare[4, 5:6]; xlabel = "[DOP] (μM)", ylabel = "z (m)", yaxisposition = :left)
# xlims!(down_DOP, -0.01, 0.1)
# lines!(down_DOP, init_down_DOPₙ, zw[1:200],color = :red, linestyle = :dash, label = "Initial")
# DOP_down = lines!(down_DOP, down_DOPₙ[][:, 1], zw[1:200], label = "Dynamic")
# axislegend(down_DOP, position = :rb)

fig_compare[1, 1:4] = Label(fig_compare, title, tellwidth=false)

# And, finally, we record a movie.
frames = 1:length(times)
record(fig_compare, "AMOC4_test.mp4", frames, framerate=50) do i
    n[] = i
    PO4_up[1] = up_PO4ₙ[][:, 1]
    POP_up[1] = up_POPₙ[][:, 1]
    # DOP_up[1] = up_DOPₙ[][:, 1]
    PO4_down[1] = down_PO4ₙ[][:, 1]
    POP_down[1] = down_POPₙ[][:, 1]
    # DOP_down[1] = down_DOPₙ[][:, 1]
end
nothing #hide

################################# PO4 & POP time-conc plot #################################
Nz = 100
time_after = 180

PO4_t_up = zeros(Nz,time_after) # zw,n
PO4_t_up_diff = zeros(Nz,time_after) 
POP_t_up = zeros(Nz,time_after) # zw,n
POP_t_up_diff = zeros(Nz,time_after) 

PO4_t_down = zeros(Nz,time_after) # zw,n
PO4_t_down_diff = zeros(Nz,time_after) 
POP_t_down = zeros(Nz,time_after) # zw,n
POP_t_down_diff = zeros(Nz,time_after) 

# add downwell
for i = 1:time_after # n
    PO4_t_up[:, i] = 1e3*interior(PO4_timeseries[i], 100, 1, :)
    PO4_t_up_diff[:, i] = 1e3*(interior(PO4_timeseries[i], 100, 1, :).-interior(PO4_timeseries[1], 100, 1, :))
    POP_t_up[:, i] = 1e3*interior(POP_timeseries[i], 100, 1, :)
    POP_t_up_diff[:, i] = 1e3*(interior(POP_timeseries[i], 100, 1, :).-interior(POP_timeseries[1], 100, 1, :))
    PO4_t_down[:, i] = 1e3*interior(PO4_timeseries[i], 440, 1, :)
    PO4_t_down_diff[:, i] = 1e3*(interior(PO4_timeseries[i], 440, 1, :).-interior(PO4_timeseries[1], 440, 1, :))
    POP_t_down[:, i] = 1e3*interior(POP_timeseries[i], 440, 1, :)
    POP_t_down_diff[:, i] = 1e3*(interior(POP_timeseries[i], 440, 1, :).-interior(POP_timeseries[1], 440, 1, :))
end

fig_time = Figure(size=(1600, 800))

ax_PO4_time = Axis(fig_time[1, 1]; xlabel = "time (days)", ylabel = "z (m)", title = "Surface [PO₄] (μM) vs time", aspect = 1)
hm_PO4_time = heatmap!(ax_PO4_time, 1:time_after, zw[85:100], PO4_t_up[85:100,:]'; colorrange = (0,1),colormap = :jet1) 
Colorbar(fig_time[1, 2], hm_PO4_time; flipaxis = false)

ax_PO4_time = Axis(fig_time[1, 3]; xlabel = "time (days)", ylabel = "z (m)", title = "Surface Δ[PO₄] (μM) vs time", aspect = 1)
hm_PO4_time = heatmap!(ax_PO4_time, 1:time_after, zw[85:100], PO4_t_up_diff[85:100,:]'; colorrange = (-0.02,0.02),colormap = :balance) 
Colorbar(fig_time[1, 4], hm_PO4_time; flipaxis = false)

ax_POP_time = Axis(fig_time[1, 5]; xlabel = "time (days)", ylabel = "z (m)", title = "[POP] (μM) vs time", aspect = 1)
hm_POP_time = heatmap!(ax_POP_time, 1:time_after, zw, POP_t_up'; colorrange = (0,0.001),colormap = :jet1) 
Colorbar(fig_time[1, 6], hm_POP_time; flipaxis = false)

ax_POP_time = Axis(fig_time[1, 7]; xlabel = "time (days)", ylabel = "z (m)", title = "Δ[POP] (μM) vs time", aspect = 1)
hm_POP_time = heatmap!(ax_POP_time, 1:time_after, zw, POP_t_up_diff'; colorrange = (-0.001,0.001),colormap = :balance) 
Colorbar(fig_time[1, 8], hm_POP_time; flipaxis = false)

ax_PO4_time2 = Axis(fig_time[2, 1]; xlabel = "time (days)", ylabel = "z (m)", title = "Surface [PO₄] (μM) vs time", aspect = 1)
hm_PO4_time2 = heatmap!(ax_PO4_time2, 1:time_after, zw[85:100], PO4_t_down[85:100,:]'; colorrange = (0,0.5),colormap = :jet1) 
Colorbar(fig_time[2, 2], hm_PO4_time2; flipaxis = false)

ax_PO4_time2 = Axis(fig_time[2, 3]; xlabel = "time (days)", ylabel = "z (m)", title = "Surface Δ[PO₄] (μM) vs time", aspect = 1)
hm_PO4_time2 = heatmap!(ax_PO4_time2, 1:time_after, zw[85:100], PO4_t_down_diff[85:100,:]'; colorrange = (-0.02,0.02),colormap = :balance) 
Colorbar(fig_time[2, 4], hm_PO4_time2; flipaxis = false)

ax_POP_time2 = Axis(fig_time[2, 5]; xlabel = "time (days)", ylabel = "z (m)", title = "[POP] (μM) vs time", aspect = 1)
hm_POP_time2 = heatmap!(ax_POP_time2, 1:time_after, zw, POP_t_down'; colorrange = (0,0.001),colormap = :jet1) 
Colorbar(fig_time[2, 6], hm_POP_time2; flipaxis = false)

ax_POP_time2 = Axis(fig_time[2, 7]; xlabel = "time (days)", ylabel = "z (m)", title = "Δ[POP] (μM) vs time", aspect = 1)
hm_POP_time2 = heatmap!(ax_POP_time2, 1:time_after, zw, POP_t_down_diff'; colorrange = (-0.0005,0.0005),colormap = :balance) 
Colorbar(fig_time[2, 8], hm_POP_time2; flipaxis = false)

display(fig_time)

################################# Plot fluxes ################################# 

# Peak flux: POP_t_up[185,] 
F100 = zeros(1,time_after) 
F200 = zeros(1,time_after) 
F500 = zeros(1,time_after) 
F1000 = zeros(1,time_after) 
F2000 = zeros(1,time_after) 

for i in 1:time_after
    F100[1,i] = 10 .* POP_t_up[95,i]
    F200[1,i] = 10 .* POP_t_up[90,i] #/ POP_t_up[93,i] 
    F500[1,i] = 10 .* POP_t_up[75,i] 
    F1000[1,i] = 10 .* POP_t_up[50,i] 
    F2000[1,i] = 10 .* POP_t_up[1,i] 
end
# z₀ = -log(0.01)*25 
# F200_martin = ones(1,time_after) .* ((205+z₀)/(150+z₀))^-0.84
# F1000_martin = ones(1,time_after) .* ((1005+z₀)/(150+z₀))^-0.84
# F2000_martin = ones(1,time_after) .* ((1995+z₀)/(150+z₀))^-0.84

fig_flux = Figure(size=(600, 600))

ax_F1000 = Axis(fig_flux[1, 1]; xlabel = "time (days)", ylabel = "F(z) (mmol m⁻² d⁻¹)", yaxisposition = :left)
# ylims!(ax_F1000, 0, 1)
lines!(ax_F1000, 1:time_after, F100[1,:], color = :orange, label = "100 m")
lines!(ax_F1000, 1:time_after, F200[1,:], color = :black, label = "200 m")
# lines!(ax_F1000, 1:time_after, F200_martin[1,:], color = :black,  linestyle = :dash, label = "200 m (Martin)")
lines!(ax_F1000, 1:time_after, F500[1,:], color = :green, label = "500 m")
lines!(ax_F1000, 1:time_after, F1000[1,:], color = :red, label = "1000 m")
# lines!(ax_F1000, 1:time_after, F1000_martin[1,:], color = :red,  linestyle = :dash, label = "1000 m (Martin)")
lines!(ax_F1000, 1:time_after, F2000[1,:], color = :blue, label = "2000 m")
# lines!(ax_F1000, 1:time_after, F2000_martin[1,:], color = :blue,  linestyle = :dash, label = "2000 m (Martin)")

axislegend(ax_F1000, position = :rt)
display(fig_flux)

########################### Plot prod/remin rate changes ############################
PO4_final = 1e3*interior(PO4_timeseries[end], :, 1, :)
POP_final = 1e3*interior(POP_timeseries[end], :, 1, :)
DOP_final = 1e3*interior(DOP_timeseries[end], :, 1, :)
Fe_final = 1e6*interior(Fe_timeseries[end], :, 1, :)
NO3_final = 1e3*interior(NO3_timeseries[end], :, 1, :)

PO4_mid= 1e3*interior(PO4_timeseries[30], :, 1, :)
POP_mid= 1e3*interior(POP_timeseries[30], :, 1, :)
DOP_mid= 1e3*interior(DOP_timeseries[30], :, 1, :)
Fe_mid= 1e6*interior(Fe_timeseries[30], :, 1, :)
NO3_mid= 1e3*interior(NO3_timeseries[30], :, 1, :)

PO4_init= 1e3*interior(PO4_timeseries[1], :, 1, :)
POP_init= 1e3*interior(POP_timeseries[1], :, 1, :)
DOP_init= 1e3*interior(DOP_timeseries[1], :, 1, :)
Fe_init= 1e6*interior(Fe_timeseries[1], :, 1, :)
NO3_init= 1e3*interior(NO3_timeseries[1], :, 1, :)

zi = zw[1:100]
z_matrix = repeat(zw[1:100], 1, 500)
μᵖ= 1e-5 # /day
kᴵ=10
kᴾ=1e-7*1024.5 * 1e3 # unit conversion 1e3
kᴺ=1.6e-6*1024.5 * 1e3 # unit conversion 1e3
kᶠ=1e-10*1024.5 * 1e6 # unit conversion 1e6
PAR_init = 700 
PAR_final = 1400 

I_begin = max.(0, PAR_init .* exp.(z_matrix' ./ 25))
I_end = max.(0, PAR_final .* exp.(z_matrix' ./ 25))
P_begin = max.(0, PO4_init)
P_mid = max.(0, PO4_mid)
P_end = max.(0, PO4_final)
N_begin = max.(0, NO3_init)
N_mid = max.(0, NO3_mid)
N_end = max.(0, NO3_final)
F_begin = max.(0, Fe_init)
F_mid = max.(0, Fe_mid)
F_end = max.(0, Fe_final)

light_lim_begin = I_begin ./ (I_begin .+ kᴵ)
light_lim_end = I_end ./ (I_end .+ kᴵ)
p_lim_begin = (P_begin ./ (P_begin .+ kᴾ))
n_lim_begin = (N_begin ./ (N_begin .+ kᴺ))
f_lim_begin = (F_begin ./ (F_begin .+ kᶠ))
p_lim_mid = (P_mid ./ (P_mid .+ kᴾ))
n_lim_mid = (N_mid ./ (N_mid .+ kᴺ))
f_lim_mid = (F_mid ./ (F_mid .+ kᶠ))
p_lim_end = (P_end ./ (P_end .+ kᴾ))
n_lim_end = (N_end ./ (N_end .+ kᴺ))
f_lim_end = (F_end ./ (F_end .+ kᶠ))

NCP_init = max.(0, μᵖ .* light_lim_begin .* min.(p_lim_begin, n_lim_begin, f_lim_begin))
NCP_mid = max.(0, μᵖ .* light_lim_end .* min.(p_lim_mid, n_lim_mid, f_lim_mid))
NCP_final = max.(0, μᵖ .* light_lim_end .* min.(p_lim_end, n_lim_end, f_lim_end))
NCP_diff = NCP_final .- NCP_init
   
z₀ = log(0.01)*25 
rₛₑ = 1 #/day
Remin_init = -0.84 .* 10 ./(z_matrix' .+ z₀) .* POP_init
Remin_init[:,1] = rₛₑ .* POP_init[:,1] # sedimentary remin

Remin_mid = -0.84 .* 10 ./(z_matrix' .+ z₀) .* POP_mid
Remin_mid[:,1] = rₛₑ .* POP_mid[:,1] # sedimentary remin

Remin_final = -0.84 .* 10 ./(z_matrix' .+ z₀) .* POP_final
# Remin_final = ifelse(z < -1990, rₛₑ .* POP_final, 0.1 .* POP_final)
Remin_final[:,1] = rₛₑ .* POP_final[:,1] # sedimentary remin

Remin_diff = Remin_final .- Remin_init

fig_rate = Figure(size = (1000, 1000))

ax_NCP = Axis(fig_rate[1, 1]; xlabel = "x (m)", ylabel = "z (m)",title = "Top 500 m ΔNCP (mmol m⁻³ d⁻¹)", aspect = 1)
hm_NCP = heatmap!(ax_NCP, xw, zw[75:100], 1000*NCP_diff[:,75:100]; colorrange = (-0.0005,0.0005),colormap = :balance) 
Colorbar(fig_rate[1, 2], hm_NCP; flipaxis = false)

ax_remin = Axis(fig_rate[1, 3]; xlabel = "x (m)", ylabel = "z (m)",title = "ΔRemin (mmol m⁻³ d⁻¹)", aspect = 1)
hm_remin = heatmap!(ax_remin, xw, zw, Remin_diff; colorrange = (-0.0000075,0.0000075),colormap = :balance) 
Colorbar(fig_rate[1, 4], hm_remin; flipaxis = false)

avg_NCP_init = mean(1000*NCP_init; dims = 1)
avg_NCP_mid = mean(1000*NCP_mid; dims = 1)
avg_NCP_final = mean(1000*NCP_final; dims = 1)

ax_avg_NCP = Axis(fig_rate[2, 1:2]; xlabel = "Rate (mmol m⁻³ d⁻¹)", ylabel = "z (m)", title = "Ave. NCP", yaxisposition = :right)
xlims!(ax_avg_NCP, -0.0001, 0.00075)
ylims!(ax_avg_NCP, -500, 0)
lines!(ax_avg_NCP, vec(avg_NCP_init),zi, color = :blue, label = "Initial")
lines!(ax_avg_NCP, vec(avg_NCP_mid),zi, color = :black, label = "After 1 month")
lines!(ax_avg_NCP, vec(avg_NCP_final),zi, color = :red, label = "After 6 months")
axislegend(ax_avg_NCP, position = :rb)

avg_R_init = mean(Remin_init; dims = 1)
avg_R_mid = mean(Remin_mid; dims = 1)
avg_R_final = mean(Remin_final; dims = 1)

ax_avg_R = Axis(fig_rate[2, 3:4]; xlabel = "Rate (mmol m⁻³ d⁻¹)", ylabel = "z (m)", title = "Ave. POP Remin", yaxisposition = :right)
xlims!(ax_avg_R, -0.000001, 0.00002)
lines!(ax_avg_R, vec(avg_R_init),zi, color = :blue, label = "Initial")
lines!(ax_avg_R, vec(avg_R_mid),zi, color = :black, label = "After 1 month")
lines!(ax_avg_R, vec(avg_R_final),zi, color = :red, label = "After 6 months")
axislegend(ax_avg_R, position = :rb)

# compare to Martin curve
# Remin_flux = vec(mean(10 .* POP_final; dims = 1))
# Martin_flux = Remin_flux[183]*((zi[183]+z₀)./(zi.+z₀)).^0.84

# ax_flux = Axis(fig_rate[3, 3:4]; xlabel = "POP flux (mmol m⁻² d⁻¹)", ylabel = "z (m)", title = "Flux comparison", yaxisposition = :right)
# xlims!(ax_flux, 0, 0.2)
# lines!(ax_flux, Remin_flux[2:200], zi[2:200], label = "Model")
# lines!(ax_flux, Martin_flux, zi, label = "Martin curve")
# axislegend(ax_flux, position = :rb)

display(fig_rate)

########################### Plot NPP, e-ratio vs time ############################

Nz = 100
time_after = 180

PO4_t = zeros(Nz,time_after) # zw,n
NO3_t = zeros(Nz,time_after) 
Fe_t = zeros(Nz,time_after) 
POP_t = zeros(Nz,time_after) 
DOP_t = zeros(Nz,time_after) 

for i = 1:time_after # n
    PO4_t[:, i] = 1e3*interior(PO4_timeseries[i], 100, 1, :)
    NO3_t[:, i] = 1e3*(interior(NO3_timeseries[i], 100, 1, :))
    Fe_t[:, i] = 1e6*(interior(Fe_timeseries[i], 100, 1, :))
    POP_t[:, i] = 1e3*interior(POP_timeseries[i], 100, 1, :)
    DOP_t[:, i] = 1e3*interior(DOP_timeseries[i], 100, 1, :)
end

zi = zw[1:100]
μᵖ= 1e-5 # /day
kᴵ=10
kᴾ=1e-7*1024.5 * 1e3 # unit conversion 1e3
kᴺ=1.6e-6*1024.5 * 1e3 # unit conversion 1e3
kᶠ=1e-10*1024.5 * 1e6 # unit conversion 1e6
PAR_init = 700 
PAR_final = 1400 

I_begin = max.(0, PAR_init .* exp.(zi ./ 25))
I_end = max.(0, PAR_final .* exp.(zi ./ 25))

I_t = zeros(Nz,time_after)
I_t[:, 1] = I_begin 
for i = 2:time_after # n
    I_t[:, i] = I_end
end
P_t = max.(0, PO4_t)
N_t = max.(0, NO3_t)
F_t = max.(0, Fe_t)

light_lim_t = I_t ./ (I_t .+ kᴵ)
p_lim_t = (P_t ./ (P_t .+ kᴾ))
n_lim_t = (N_t ./ (N_t .+ kᴺ))
f_lim_t = (F_t ./ (F_t .+ kᶠ))

NCP_t = max.(0, μᵖ .* light_lim_t .* min.(p_lim_t, n_lim_t, f_lim_t))
   
z₀ = log(0.01)*25 
rₛₑ = 1 #/day
z_matrix = repeat(zi, 1, time_after)
Remin_t = -0.84 .* 10 ./(z_matrix .+ z₀) .* POP_t
Remin_t[1,:] = rₛₑ .* POP_t[1,:] # sedimentary remin
DOPRemin = 1/30 .* DOP_t

POP_flux = 10 .* POP_t # mmol m-2 d-1

# Calculate e-ratio (250 m, zw[88])
NCP_sum = zeros(1,180)
topR_sum = zeros(1,180)
upperR_sum = zeros(1,180)
btmR_sum = zeros(1,180)
# DOPR_sum = zeros(1,180)
for t_index = 1:180
    NCP_sum[t_index] = 1e3 .* sum(NCP_t[:,t_index])
    topR_sum[t_index] = sum(Remin_t[1:88,t_index]) # Below 250 m 
    upperR_sum[t_index] = sum(Remin_t[1:50,t_index]) # Below 1000 m 
end
eratio = topR_sum./ NCP_sum
SeqEff = upperR_sum./ NCP_sum

fig_PR = Figure(size = (600, 1000))

ax_NCP = Axis(fig_PR[1, 1]; xlabel = "time (days)", ylabel = "z (m)",title = "Top 500 m NCP (mmol m⁻³ d⁻¹)", aspect = 1)
hm_NCP = heatmap!(ax_NCP, 1:time_after, zw[75:100], 1000*NCP_t[75:100,:]'; colorrange = (0,1.5e-3),colormap = :jet1) 
Colorbar(fig_PR[1, 2], hm_NCP; flipaxis = false)

ax_remin = Axis(fig_PR[2, 1]; xlabel = "time (days)", ylabel = "z (m)",title = "POP remineralization (mmol m⁻³ d⁻¹)", aspect = 1)
hm_remin = heatmap!(ax_remin, 1:time_after, zw, Remin_t'; colorrange = (0,3e-5),colormap = :jet1) 
Colorbar(fig_PR[2, 2], hm_remin; flipaxis = false)

ax_eRatio = Axis(fig_PR[3, 1:2]; xlabel = "time (days)", ylabel = "Ratio", title = "Effeciency", yaxisposition = :left)
ylims!(ax_eRatio, 0, 0.5)
lines!(ax_eRatio, 1:time_after, vec(eratio), color = :black, label = "Export effeciency (250 m)")
lines!(ax_eRatio, 1:time_after, vec(SeqEff), color = :blue, label = "Sequestration effeciency (1000 m)")
axislegend(ax_eRatio, position = :rt)

display(fig_PR)
