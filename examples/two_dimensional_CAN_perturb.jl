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
Nx = 100 
Nz = 200
Lx = 1000kilometers   # m
Lz = 2000           # m

# We use a two-dimensional grid, with a `Flat` `y`-direction:
grid = RectilinearGrid(#GPU(),
                       size = (Nx, Nz),
                       x = (0, Lx),
                       z = (-Lz, 0),
                       topology=(Bounded, Flat, Bounded))

#################################### Boundary conditions ####################################

const τ₀ = 2e-4           # m² s⁻²
function τy(x, t, p)
    xfrac = x / p.Lx
    return xfrac < 0.8 ?  p.τ₀ * sinpi(xfrac / 1.6) :  p.τ₀ * sinpi((1 - xfrac) / 0.4)
end
y_wind_stress = FluxBoundaryCondition(τy, parameters=(; τ₀=τ₀, Lx=Lx)) 
v_bcs = FieldBoundaryConditions(top=y_wind_stress)

#################################### Model ####################################

# mixing_length = CATKEMixingLength(Cᵇ = 0.001)
# catke = CATKEVerticalDiffusivity(; mixing_length, tke_time_step = 10minutes)

kz(x,z,t) = 1e-5 + 5e-4 * (tanh((z+75)/10)+1) + 1e-3 * exp(-(z+2000)/50)
vertical_closure = VerticalScalarDiffusivity(;ν=kz, κ=kz)
horizontal_closure = HorizontalScalarDiffusivity(ν=1e3, κ=1e3)

M² = 0         # 1e-7 s⁻², squared buoyancy frequency
N² = 0         # 1e-5 s⁻²
bᵢ(x, z) = N² * z + M² * x

############################ Reset nutrient concentrations ############################ 

using JLD2
data = JLD2.load("CAN_2D_Rz17.jld2")
pop_keys = filter(key -> startswith(key, "timeseries/POP/") && occursin(r"\d+$", key), keys(data))
numbers = [parse(Int, match(r"\d+$", key).match) for key in pop_keys]
max_number = maximum(numbers)

mid_POP = data["timeseries/POP/$max_number"]
mid_DOP = data["timeseries/DOP/$max_number"]
mid_PO₄ = data["timeseries/PO₄/$max_number"]
mid_NO₃ = data["timeseries/NO₃/$max_number"]
mid_Fe = data["timeseries/Fe/$max_number"]
mid_DIC = data["timeseries/DIC/$max_number"]
mid_ALK = data["timeseries/ALK/$max_number"]

mid_u = data["timeseries/u/$max_number"]
mid_v = data["timeseries/v/$max_number"]
mid_w = data["timeseries/w/$max_number"]

# Set deep nutrient concentrations 
# mid_PO₄[95:100, 1, 195:200] .= 5e-3

model2 = HydrostaticFreeSurfaceModel(; grid,
                                    buoyancy = BuoyancyTracer(),
                                    biogeochemistry = CarbonAlkalinityNutrients(; grid,
                                                                        maximum_net_community_production_rate  = 1e-4/day,
                                                                        incident_PAR  = 1400,
                                                                        iron_scavenging_rate = 0),
                                    coriolis = FPlane(; f=1e-4),
                                    closure = (vertical_closure, horizontal_closure), 
                                    tracers = (:b, :e, :DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe),
                                    momentum_advection = WENO(),
                                    tracer_advection = WENO(),
                                    boundary_conditions = (; v=v_bcs)) 

set!(model2, b=bᵢ, DIC=mid_DIC, ALK=mid_ALK, NO₃=mid_NO₃, PO₄=mid_PO₄, DOP=mid_DOP, POP=mid_POP, Fe = mid_Fe) 

simulation2 = Simulation(model2; Δt = 5minutes, stop_time=365.25days) 
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
                     filename = "CAN_2D_Rz17_perturbPAR",
                     overwrite_existing = true)

run!(simulation2)

############################ Visualizing the solution ############################

filepath = simulation2.output_writers[:simple_output].filepath
# filepath = "./CAN_2D_Rz12_reset.jld2"

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
hm_u = heatmap!(ax_u, xw, zw, uₙ; colorrange = (-0.5, 0.5), colormap = :balance) 
Colorbar(fig[2, 2], hm_u; flipaxis = false)

ax_w = Axis(fig[2, 3]; xlabel = "x (m)", ylabel = "z (m)", title = "w (cm/s)", aspect = 1)
hm_w = heatmap!(ax_w, xw, zw, wₙ; colorrange = (-2e-3,2e-3), colormap = :balance) 
Colorbar(fig[2, 4], hm_w; flipaxis = false)

ax_PO4 = Axis(fig[3, 1]; xlabel = "x (m)", ylabel = "z (m)", title = "[PO₄] (μM)", aspect = 1)
hm_PO4 = heatmap!(ax_PO4, xw, zw, PO4ₙ; colormap = :jet1) 
Colorbar(fig[3, 2], hm_PO4; flipaxis = false)

ax_POP = Axis(fig[3, 3]; xlabel = "x (m)", ylabel = "z (m)", title = "[POP] (μM)", aspect = 1)
hm_POP = heatmap!(ax_POP, xw, zw, POPₙ; colorrange = (0,0.1),colormap = :jet1) 
Colorbar(fig[3, 4], hm_POP; flipaxis = false)

ax_avg_PO4 = Axis(fig[4, 1:2]; xlabel = "[PO₄] (μM)", ylabel = "z (m)", title = "Average [PO₄] (μM)", yaxisposition = :right)
xlims!(ax_avg_PO4, 0, 7)
PO4_prof = lines!(ax_avg_PO4, avg_PO4ₙ[][1, :], zw[1:200])

ax_avg_POP = Axis(fig[4, 3:4]; xlabel = "[POP] (μM)", ylabel = "z (m)", title = "Average [POP] (μM)",yaxisposition = :right)
xlims!(ax_avg_POP, 0, 0.04)
POP_prof = lines!(ax_avg_POP, avg_POPₙ[][1, :], zw[1:200])

fig[1, 1:4] = Label(fig, title, tellwidth=false)

# And, finally, we record a movie.
frames = 1:length(times)
record(fig, "CAN_2D_Rz17_perturbPAR.mp4", frames, framerate=50) do i
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
DOP_init = 1e3*interior(DOP_timeseries[1], :, 1, :)
# sum(PO4_init.+POP_init.+DOP_init)

n = Observable(1)
title = @lift @sprintf("t = Day %d", times[$n] / 1days) 

diff_PO4 = @lift 1e3*(interior(PO4_timeseries[$n], :, 1, :) .- interior(PO4_timeseries[1], :, 1, :))
diff_POP = @lift 1e3*(interior(POP_timeseries[$n], :, 1, :) .- interior(POP_timeseries[1], :, 1, :))
diff_DOP = @lift 1e3*(interior(DOP_timeseries[$n], :, 1, :) .- interior(DOP_timeseries[1], :, 1, :))
# sum(diff_PO4)

# One upwelling station
up_PO4ₙ = @lift 1e3*interior(PO4_timeseries[$n], 97, 1, :)
up_DOPₙ = @lift 1e3*interior(DOP_timeseries[$n], 97, 1, :)
up_POPₙ = @lift 1e3*interior(POP_timeseries[$n], 97, 1, :)

init_up_PO4ₙ = 1e3*interior(PO4_timeseries[1], 97, 1, :)
init_up_POPₙ = 1e3*interior(POP_timeseries[1], 97, 1, :)
init_up_DOPₙ = 1e3*interior(DOP_timeseries[1], 97, 1, :)

# One downwelling station
down_PO4ₙ = @lift 1e3*interior(PO4_timeseries[$n], 40, 1, :)
down_DOPₙ = @lift 1e3*interior(DOP_timeseries[$n], 40, 1, :)
down_POPₙ = @lift 1e3*interior(POP_timeseries[$n], 40, 1, :)

init_down_PO4ₙ = 1e3*interior(PO4_timeseries[1], 40, 1, :)
init_down_POPₙ = 1e3*interior(POP_timeseries[1], 40, 1, :)
init_down_DOPₙ = 1e3*interior(DOP_timeseries[1], 40, 1, :)

fig_compare = Figure(size=(1100, 1100))

# Row 1: final - initial
ax_PO4 = Axis(fig_compare[2, 1]; xlabel = "x (m)", ylabel = "z (m)", title = "Δ[PO₄] (μM)", aspect = 1)
hm_PO4 = heatmap!(ax_PO4, xw, zw, diff_PO4; colorrange = (-0.2,0.2),colormap = :balance) 
Colorbar(fig_compare[2, 2], hm_PO4; flipaxis = false)

ax_POP = Axis(fig_compare[2, 3]; xlabel = "x (m)", ylabel = "z (m)", title = "Δ[POP] (μM)", aspect = 1)
hm_POP = heatmap!(ax_POP, xw, zw, diff_POP; colorrange = (-0.01,0.01),colormap = :balance) 
Colorbar(fig_compare[2, 4], hm_POP; flipaxis = false)

ax_DOP = Axis(fig_compare[2, 5]; xlabel = "x (m)", ylabel = "z (m)", title = "Δ[DOP] (μM)", aspect = 1)
hm_DOP = heatmap!(ax_DOP, xw, zw, diff_DOP; colorrange = (-0.06,0.06),colormap = :balance) 
Colorbar(fig_compare[2, 6], hm_DOP; flipaxis = false)

# Row 2: upwelling site
up_PO4 = Axis(fig_compare[3, 1:2]; xlabel = "[PO₄] (μM)", ylabel = "z (m)", yaxisposition = :left)
xlims!(up_PO4, 0, 7)
lines!(up_PO4, init_up_PO4ₙ, zw[1:200],color = :red, linestyle = :dash, label = "Initial")
PO4_up = lines!(up_PO4, up_PO4ₙ[][:, 1], zw[1:200], label = "Dynamic")
axislegend(up_PO4, position = :lb)

up_POP = Axis(fig_compare[3, 3:4]; xlabel = "[POP] (μM)", ylabel = "z (m)", yaxisposition = :left)
xlims!(up_POP, 0, 0.1)
lines!(up_POP, init_up_POPₙ, zw[1:200],color = :red, linestyle = :dash, label = "Initial")
POP_up = lines!(up_POP, up_POPₙ[][:, 1], zw[1:200], label = "Dynamic")
axislegend(up_POP, position = :rb)

up_DOP = Axis(fig_compare[3, 5:6]; xlabel = "[DOP] (μM)", ylabel = "z (m)", yaxisposition = :left)
xlims!(up_DOP, -0.02, 0.5)
lines!(up_DOP, init_up_DOPₙ, zw[1:200],color = :red, linestyle = :dash, label = "Initial")
DOP_up = lines!(up_DOP, up_DOPₙ[][:, 1], zw[1:200], label = "Dynamic")
axislegend(up_DOP, position = :rb)

# Row 3: downwelling site
down_PO4 = Axis(fig_compare[4, 1:2]; xlabel = "[PO₄] (μM)", ylabel = "z (m)", yaxisposition = :left)
xlims!(down_PO4, 0, 7)
lines!(down_PO4, init_down_PO4ₙ, zw[1:200],color = :red, linestyle = :dash, label = "Initial")
PO4_down = lines!(down_PO4, down_PO4ₙ[][:, 1], zw[1:200], label = "Dynamic")
axislegend(down_PO4, position = :lb)

down_POP = Axis(fig_compare[4, 3:4]; xlabel = "[POP] (μM)", ylabel = "z (m)", yaxisposition = :left)
xlims!(down_POP, 0, 0.01)
lines!(down_POP, init_down_POPₙ, zw[1:200],color = :red, linestyle = :dash, label = "Initial")
POP_down = lines!(down_POP, down_POPₙ[][:, 1], zw[1:200], label = "Dynamic")
axislegend(down_POP, position = :rb)

down_DOP = Axis(fig_compare[4, 5:6]; xlabel = "[DOP] (μM)", ylabel = "z (m)", yaxisposition = :left)
xlims!(down_DOP, -0.01, 0.1)
lines!(down_DOP, init_down_DOPₙ, zw[1:200],color = :red, linestyle = :dash, label = "Initial")
DOP_down = lines!(down_DOP, down_DOPₙ[][:, 1], zw[1:200], label = "Dynamic")
axislegend(down_DOP, position = :rb)

fig_compare[1, 1:6] = Label(fig_compare, title, tellwidth=false)

# And, finally, we record a movie.
frames = 1:length(times)
record(fig_compare, "CAN_2D_Rz17_perturbPAR_compare.mp4", frames, framerate=75) do i
    n[] = i
    PO4_up[1] = up_PO4ₙ[][:, 1]
    POP_up[1] = up_POPₙ[][:, 1]
    DOP_up[1] = up_DOPₙ[][:, 1]
    PO4_down[1] = down_PO4ₙ[][:, 1]
    POP_down[1] = down_POPₙ[][:, 1]
    DOP_down[1] = down_DOPₙ[][:, 1]
end
nothing #hide

################################# PO4 & POP time-conc plot #################################
PO4_t_up = zeros(200,366) # zw,n
PO4_t_up_diff = zeros(200,366) 
POP_t_up = zeros(200,366) # zw,n
POP_t_up_diff = zeros(200,366) 

PO4_t_down = zeros(200,366) # zw,n
PO4_t_down_diff = zeros(200,366) 
POP_t_down = zeros(200,366) # zw,n
POP_t_down_diff = zeros(200,366) 

# add downwell
for i = 1:366 # n
    PO4_t_up[:, i] = 1e3*interior(PO4_timeseries[i], 94, 1, :)
    PO4_t_up_diff[:, i] = 1e3*(interior(PO4_timeseries[i], 94, 1, :).-interior(PO4_timeseries[1], 94, 1, :))
    POP_t_up[:, i] = 1e3*interior(POP_timeseries[i], 94, 1, :)
    POP_t_up_diff[:, i] = 1e3*(interior(POP_timeseries[i], 94, 1, :).-interior(POP_timeseries[1], 94, 1, :))
    PO4_t_down[:, i] = 1e3*interior(PO4_timeseries[i], 40, 1, :)
    PO4_t_down_diff[:, i] = 1e3*(interior(PO4_timeseries[i], 40, 1, :).-interior(PO4_timeseries[1], 40, 1, :))
    POP_t_down[:, i] = 1e3*interior(POP_timeseries[i], 40, 1, :)
    POP_t_down_diff[:, i] = 1e3*(interior(POP_timeseries[i], 40, 1, :).-interior(POP_timeseries[1], 40, 1, :))
end

fig_time = Figure(size=(1600, 800))

ax_PO4_time = Axis(fig_time[1, 1]; xlabel = "time (days)", ylabel = "z (m)", title = "Surface [PO₄] (μM) vs time", aspect = 1)
hm_PO4_time = heatmap!(ax_PO4_time, 1:366, zw[185:200], PO4_t_up[185:200,:]'; colorrange = (0,4),colormap = :jet1) 
Colorbar(fig_time[1, 2], hm_PO4_time; flipaxis = false)

ax_PO4_time = Axis(fig_time[1, 3]; xlabel = "time (days)", ylabel = "z (m)", title = "Surface Δ[PO₄] (μM) vs time", aspect = 1)
hm_PO4_time = heatmap!(ax_PO4_time, 1:366, zw[185:200], PO4_t_up_diff[185:200,:]'; colorrange = (-3,3),colormap = :balance) 
Colorbar(fig_time[1, 4], hm_PO4_time; flipaxis = false)

ax_POP_time = Axis(fig_time[1, 5]; xlabel = "time (days)", ylabel = "z (m)", title = "[POP] (μM) vs time", aspect = 1)
hm_POP_time = heatmap!(ax_POP_time, 1:366, zw, POP_t_up'; colorrange = (0,0.12),colormap = :jet1) 
Colorbar(fig_time[1, 6], hm_POP_time; flipaxis = false)

ax_POP_time = Axis(fig_time[1, 7]; xlabel = "time (days)", ylabel = "z (m)", title = "Δ[POP] (μM) vs time", aspect = 1)
hm_POP_time = heatmap!(ax_POP_time, 1:366, zw, POP_t_up_diff'; colorrange = (-0.05,0.05),colormap = :balance) 
Colorbar(fig_time[1, 8], hm_POP_time; flipaxis = false)

ax_PO4_time2 = Axis(fig_time[2, 1]; xlabel = "time (days)", ylabel = "z (m)", title = "Surface [PO₄] (μM) vs time", aspect = 1)
hm_PO4_time2 = heatmap!(ax_PO4_time2, 1:366, zw[185:200], PO4_t_down[185:200,:]'; colorrange = (0,1),colormap = :jet1) 
Colorbar(fig_time[2, 2], hm_PO4_time2; flipaxis = false)

ax_PO4_time2 = Axis(fig_time[2, 3]; xlabel = "time (days)", ylabel = "z (m)", title = "Surface Δ[PO₄] (μM) vs time", aspect = 1)
hm_PO4_time2 = heatmap!(ax_PO4_time2, 1:366, zw[185:200], PO4_t_down_diff[185:200,:]'; colorrange = (-1,1),colormap = :balance) 
Colorbar(fig_time[2, 4], hm_PO4_time2; flipaxis = false)

ax_POP_time2 = Axis(fig_time[2, 5]; xlabel = "time (days)", ylabel = "z (m)", title = "[POP] (μM) vs time", aspect = 1)
hm_POP_time2 = heatmap!(ax_POP_time2, 1:366, zw, POP_t_down'; colorrange = (0,0.02),colormap = :jet1) 
Colorbar(fig_time[2, 6], hm_POP_time2; flipaxis = false)

ax_POP_time2 = Axis(fig_time[2, 7]; xlabel = "time (days)", ylabel = "z (m)", title = "Δ[POP] (μM) vs time", aspect = 1)
hm_POP_time2 = heatmap!(ax_POP_time2, 1:366, zw, POP_t_down_diff'; colorrange = (-0.0025,0.0025),colormap = :balance) 
Colorbar(fig_time[2, 8], hm_POP_time2; flipaxis = false)

display(fig_time)

# Flux (zw[190] = 100 m, zw[100] = 1,000 m)
fig_flux = Figure(size=(600, 600))
# Peak flux: POP_t_up[185,] 
F200 = zeros(1,366) 
F1000 = zeros(1,366) 
F2000 = zeros(1,366) 
for i in 1:366
    F200[1,i] = POP_t_up[180,i] / POP_t_up[185,i] 
    F1000[1,i] = POP_t_up[100,i] / POP_t_up[185,i] 
    F2000[1,i] = POP_t_up[1,i] / POP_t_up[185,i] 
end
z₀ = -log(0.01)*25 
F200_martin = ones(1,366) .* ((205+z₀)/(150+z₀))^-0.84
F1000_martin = ones(1,366) .* ((1005+z₀)/(150+z₀))^-0.84
F2000_martin = ones(1,366) .* ((1995+z₀)/(150+z₀))^-0.84

ax_F1000 = Axis(fig_flux[1, 1]; xlabel = "time (days)", ylabel = "F(z)/F₁₅₀ₘ", yaxisposition = :left)
ylims!(ax_F1000, 0, 1)
lines!(ax_F1000, 1:366, F200[1,:], color = :black, label = "200 m")
# lines!(ax_F1000, 1:366, F200_martin[1,:], color = :black,  linestyle = :dash, label = "200 m (Martin)")
lines!(ax_F1000, 1:366, F1000[1,:], color = :red, label = "1000 m")
# lines!(ax_F1000, 1:366, F1000_martin[1,:], color = :red,  linestyle = :dash, label = "1000 m (Martin)")
lines!(ax_F1000, 1:366, F2000[1,:], color = :blue, label = "2000 m")
# lines!(ax_F1000, 1:366, F2000_martin[1,:], color = :blue,  linestyle = :dash, label = "2000 m (Martin)")

axislegend(ax_F1000, position = :rt)
display(fig_flux)

################################## Plot all BGC tracers ##################################

# data = JLD2.load("path/to/yourfile.jld2")
# ave_POP = data["timeseries/avg_POP/1827"]
#=
PO4_final = 1e3*interior(PO4_timeseries[end], :, 1, :)
POP_final = 1e3*interior(POP_timeseries[end], :, 1, :)
DOP_final = 1e3*interior(DOP_timeseries[end], :, 1, :)
Fe_final = 1e6*interior(Fe_timeseries[end], :, 1, :)
NO3_final = 1e3*interior(NO3_timeseries[end], :, 1, :)

z_matrix = repeat(zi, 1, 100)
μᵖ= 1e-4 # /day
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
rₛₑ = 1 #/day
Remin_final = -0.84 .* 10 ./(z_matrix' .+ z₀) .* POP_final
# Remin_final = ifelse(z < -1990, rₛₑ .* POP_final, 0.1 .* POP_final)

Remin_final[:,1] = rₛₑ .* POP_final[:,1] # sedimentary remin

# tracer concentration and limitation terms
fig_can = Figure(size = (1000, 1000))

# Light limitation
ax_l_lim = Axis(fig_can[1, 1]; xlabel = "x (m)", ylabel = "z (m)",title = "I/(I+kᵢ)", aspect = 1)
hm_l_lim = heatmap!(ax_l_lim, xw, zw, light_lim; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_can[1, 2], hm_l_lim; flipaxis = false)
# PO4 concentrations
ax_PO4 = Axis(fig_can[1, 3]; xlabel = "x (m)", ylabel = "z (m)", title = "[PO₄] (μM)", aspect = 1)
hm_PO4 = heatmap!(ax_PO4, xw, zw, PO4_final; colorrange = (0, 3),colormap = :jet1) 
Colorbar(fig_can[1, 4], hm_PO4; flipaxis = false)
# PO4 limitation
ax_p_lim = Axis(fig_can[1, 5]; xlabel = "x (m)", ylabel = "z (m)",title = "PO₄/(PO₄+kₚₒ₄)", aspect = 1)
hm_p_lim = heatmap!(ax_p_lim, xw, zw, p_lim; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_can[1, 6], hm_p_lim; flipaxis = false)
# POP concentrations
ax_POP = Axis(fig_can[2, 1]; xlabel = "x (m)", ylabel = "z (m)", title = "[POP] (μM)", aspect = 1)
hm_POP = heatmap!(ax_POP, xw, zw, POP_final; colorrange = (0,0.01),colormap = :jet1) 
Colorbar(fig_can[2, 2], hm_POP; flipaxis = false)
# NO3 concentrations
ax_NO3 = Axis(fig_can[2, 3]; xlabel = "x (m)", ylabel = "z (m)", title = "[NO₃] (μM)", aspect = 1)
hm_NO3 = heatmap!(ax_NO3, xw, zw, NO3_final; colorrange = (0,50), colormap = :jet1) 
Colorbar(fig_can[2, 4], hm_NO3; flipaxis = false)
# NO3 limitation
ax_n_lim = Axis(fig_can[2, 5]; xlabel = "x (m)", ylabel = "z (m)",title = "NO₃/(NO₃+kₙₒ₃)", aspect = 1)
hm_n_lim = heatmap!(ax_n_lim, xw, zw, n_lim; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_can[2, 6], hm_n_lim; flipaxis = false)
# DOP concentrations
ax_DOP = Axis(fig_can[3, 1]; xlabel = "x (m)", ylabel = "z (m)", title = "[DOP] (μM)", aspect = 1)
hm_DOP = heatmap!(ax_DOP, xw, zw, DOP_final; colorrange = (0,0.07),colormap = :jet1) 
Colorbar(fig_can[3, 2], hm_DOP; flipaxis = false)
# Fe concentrations
ax_Fe = Axis(fig_can[3, 3]; xlabel = "x (m)", ylabel = "z (m)", title = "[Fe] (nM)", aspect = 1)
hm_Fe = heatmap!(ax_Fe, xw, zw, Fe_final; colorrange = (0,1.5), colormap = :jet1) 
Colorbar(fig_can[3, 4], hm_Fe; flipaxis = false)
# Fe limitation
ax_f_lim = Axis(fig_can[3, 5]; xlabel = "x (m)", ylabel = "z (m)",title = "Fe/(Fe+kᵢᵣₒₙ)", aspect = 1)
hm_f_lim = heatmap!(ax_f_lim, xw, zw, f_lim; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_can[3, 6], hm_f_lim; flipaxis = false)

display(fig_can)

################################## Plot prod/remin rates ##################################

fig_rate = Figure(size = (800, 1000))

ax_NCP = Axis(fig_rate[1, 1]; xlabel = "x (m)", ylabel = "z (m)",title = "NCP (mmol m⁻³ d⁻¹)", aspect = 1)
hm_NCP = heatmap!(ax_NCP, xw, zw, 1000*NCP_final; colorrange = (0,0.07),colormap = :jet1) 
Colorbar(fig_rate[1, 2], hm_NCP; flipaxis = false)

ax_remin = Axis(fig_rate[1, 3]; xlabel = "x (m)", ylabel = "z (m)",title = "Remin (mmol m⁻³ d⁻¹)", aspect = 1)
hm_remin = heatmap!(ax_remin, xw, zw, Remin_final; colorrange = (0,0.0015),colormap = :jet1) 
Colorbar(fig_rate[1, 4], hm_remin; flipaxis = false)

# ax_logNCP = Axis(fig_rate[2, 1]; xlabel = "x (m)", ylabel = "z (m)",title = "Log(NCP) (mmol m⁻³ d⁻¹)", aspect = 1)
# hm_logNCP = heatmap!(ax_logNCP, xw, zw[170:200], log.(abs.(1000*NCP_final[:,170:200])); colorrange = (-4,1),colormap = :matter) 
# Colorbar(fig_rate[2, 2], hm_logNCP; flipaxis = false)

# ax_logremin = Axis(fig_rate[2, 3]; xlabel = "x (m)", ylabel = "z (m)",title = "Log(Remin) (mmol m⁻³ d⁻¹)", aspect = 1)
# hm_logremin = heatmap!(ax_logremin, xw, zw, log.(abs.(1000*Remin_final)); colorrange = (-5,2),colormap = :matter) 
# Colorbar(fig_rate[2, 4], hm_logremin; flipaxis = false)

avg_NCP = mean(1000*NCP_final; dims = 1)
ax_avg_NCP = Axis(fig_rate[2, 1:2]; xlabel = "Rate (mmol m⁻³ d⁻¹)", ylabel = "z (m)", title = "Ave. NCP", yaxisposition = :right)
xlims!(ax_avg_NCP, -0.001, 0.04)
lines!(ax_avg_NCP, vec(avg_NCP),zi)

avg_R = mean(Remin_final; dims = 1)
ax_avg_R = Axis(fig_rate[2, 3:4]; xlabel = "Rate (mmol m⁻³ d⁻¹)", ylabel = "z (m)", title = "Ave. POP Remin", yaxisposition = :right)
xlims!(ax_avg_R, -0.0001, 0.0005)
lines!(ax_avg_R, vec(avg_R),zi)

# compare to Martin curve
Remin_flux = vec(mean(10 .* POP_final; dims = 1))
Martin_flux = Remin_flux[183]*((zi[183]+z₀)./(zi.+z₀)).^0.84

ax_flux = Axis(fig_rate[3, 3:4]; xlabel = "POP flux (mmol m⁻² d⁻¹)", ylabel = "z (m)", title = "Flux comparison", yaxisposition = :right)
xlims!(ax_flux, 0, 0.2)
lines!(ax_flux, Remin_flux[2:200], zi[2:200], label = "Model")
lines!(ax_flux, Martin_flux, zi, label = "Martin curve")
axislegend(ax_flux, position = :rb)

display(fig_rate)
=#

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
