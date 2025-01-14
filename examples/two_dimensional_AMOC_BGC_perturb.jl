# Pick up from a baseline model, run 1-year perturbation
using GLMakie
using Printf
using Statistics

using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients
using Oceananigans
using Oceananigans.Models.HydrostaticFreeSurfaceModels:
                    HydrostaticFreeSurfaceModel,
                    PrescribedVelocityFields
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization
using Oceananigans.Units
using Oceananigans.Fields: ZeroField, CenterField
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans: TendencyCallsite

const Ny = 500 
const Nz = 100
const Ly = 15000kilometers   # m
const Lz = 2000           # m

arch = CPU()
# We use a two-dimensional grid, with a `Flat` `y`-direction:
grid = RectilinearGrid(arch,
                       size = (Ny, Nz),
                       y = (0, Ly),
                       z = (-Lz, 0),
                       topology=(Flat, Bounded, Bounded))

############################## Physical circulation ############################## 
# Set streamfunction
deltaN = 200kilometers   # North downwelling width
deltaS = 500kilometers  # South upwelling width
deltaZ = -1000     # Vertical asymmetry
Ψᵢ(y, z)  = -1.2 * ((1 - exp(-y / deltaS)) * (1 - exp(-(Ly - y) / deltaN)) * 
            sinpi(z / Lz) * exp(z/deltaZ))

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

# (I have to specify u to allow CheckPointer)
u = XFaceField(grid) 
fill_halo_regions!(u, arch)

############################# Diffusivity ############################# 

kz(y,z,t) = 5e-4 + 5e-3 * (tanh((z+150)/20)+1) + 1e-2 * exp(-(z+2000)/50)
tracer_vertical_closure = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), 
                                        κ=(DIC=kz,ALK=kz,PO₄=kz,NO₃=kz,DOP=kz,POP=kz,Fe=0))
tracer_horizontal_closure = HorizontalScalarDiffusivity(
                                    κ=(DIC=1e3,ALK=1e3,PO₄=1e3,NO₃=1e3,DOP=1e3,POP=1e3,Fe=0))

# Fe field & forcing
# @inline Feᵢ(y,z) = (0.02+(-z/2000)) /1e6
# Fe_input(y, z, t) = (z > -30 && y < 100kilometers ? 1e-12 : 0.0)
# F_forcing = Forcing(Fe_input)

maximum_net_community_production_rate = CenterField(grid) 
maxNCP(y,z) = (1e-5 + sinpi(y/Ly) * 5e-5)/day
set!(maximum_net_community_production_rate, maxNCP)   
fill_halo_regions!(maximum_net_community_production_rate, arch)

############################### Set PAR(y,t) ############################### 
incident_PAR = CenterField(grid) 

function seasonal_PAR(y, z, t)
    # Constants
    day_in_year = 365.25
    axial_tilt = 23.5  # Earth's axial tilt in degrees
    latitude = (y/Ly - 0.5) * 180.0  # Convert y to latitude in degrees
    
    # Solar declination (varies over the year due to Earth's tilt)
    solar_declination = axial_tilt * sinpi(2 * t / day_in_year)
    
    # Solar angle of incidence (latitude and solar declination combined)
    solar_angle = latitude - solar_declination
    
    # PAR depends on the cosine of the solar angle 
    angle_effect = cosd(solar_angle)
    angle_effect = max(0.0, angle_effect)  # Ensure no negative values
    
    # Scale PAR by angle effect
    PAR = 700 * angle_effect
    return PAR
end
# TODO: error
set!(incident_PAR, seasonal_PAR) 
fill_halo_regions!(incident_PAR, arch)

############################## Model 2: Perturbation ############################## 
model2 = HydrostaticFreeSurfaceModel(grid = grid,
                                    biogeochemistry = CarbonAlkalinityNutrients(; grid,
                                                                        maximum_net_community_production_rate  = maximum_net_community_production_rate,
                                                                        incident_PAR = incident_PAR,
                                                                        particulate_organic_phosphorus_sedremin_timescale = 0.5 / day),
                                                                        #iron_scavenging_rate = 0),
                                    velocities = PrescribedVelocityFields(; u, v, w),
                                    tracers = (:DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe),
                                    # forcing=(; Fe=F_forcing),
                                    coriolis = nothing,
                                    buoyancy = nothing,
                                    closure = (tracer_vertical_closure, tracer_horizontal_closure))

set!(model2, DIC=2.1, ALK=2.35, NO₃=2e-2, PO₄=1.4e-3, DOP=0, POP=0, Fe = 6e-7) # mol PO₄ m⁻³

spinup_time = 365.25*5000days
perturbation_time = 365days
simulation2 = Simulation(model2; Δt = 1days, stop_time=spinup_time+perturbation_time) 

# Define a callback to zero out Fe tendency
function modify_tendency!(model2)
        model2.timestepper.Gⁿ.Fe .= 0
        # model2.timestepper.Gⁿ.Fe[1,1:50,100] .= 1e-11
        return nothing
end                                        
simulation2.callbacks[:modify_Fe] = Callback(modify_tendency!, 
                                            callsite = TendencyCallsite())

outputs = (v = model2.velocities.v,
            w = model2.velocities.w,
            PO₄= model2.tracers.PO₄,
            DOP = model2.tracers.DOP,
            POP = model2.tracers.POP,
            NO₃ = model2.tracers.NO₃,
            Fe = model2.tracers.Fe)

simulation2.output_writers[:simple_output] =
            JLD2OutputWriter(model2, outputs, 
                            schedule = TimeInterval(1days), 
                            filename = "AMOC64_perturb",
                            overwrite_existing = true)

simulation2.output_writers[:checkpointer] = Checkpointer(model2,
            schedule = TimeInterval(perturbation_time),
            prefix = "AMOC64_checkpoint",
            overwrite_existing = false)

run!(simulation2, pickup = true)

############################## Model 3: New equilibrium ############################## 
#=
model3 = HydrostaticFreeSurfaceModel(grid = grid,
                                    biogeochemistry = CarbonAlkalinityNutrients(; grid,
                                                                        maximum_net_community_production_rate  = 2e-5/day,
                                                                        # incident_PAR = incident_PAR,
                                                                        particulate_organic_phosphorus_sedremin_timescale = 0.5 / day),
                                                                        #iron_scavenging_rate = 0),
                                    velocities = PrescribedVelocityFields(; u, v, w),
                                    tracers = (:DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe),
                                    # forcing=(; Fe=F_forcing),
                                    coriolis = nothing,
                                    buoyancy = nothing,
                                    closure = (tracer_vertical_closure, tracer_horizontal_closure))

set!(model2, DIC=2.1, ALK=2.35, NO₃=2.5e-2, PO₄=1.5e-3, DOP=0, POP=0, Fe = Feᵢ) # mol PO₄ m⁻³

spinup_time = 365.25*2000days
perturbation_time = 10days
new_eq_time = 365.25days
simulation3 = Simulation(model3; Δt = 1days, stop_time=spinup_time+perturbation_time+new_eq_time) 

# Define a callback to zero out Fe tendency
function modify_tendency!(model3)
        model3.timestepper.Gⁿ.Fe .= 0
        return nothing
end                                        
simulation3.callbacks[:modify_Fe] = Callback(modify_tendency!, 
                                            callsite = TendencyCallsite())

outputs = (v = model3.velocities.v,
            w = model3.velocities.w,
            PO₄= model3.tracers.PO₄,
            DOP = model3.tracers.DOP,
            POP = model3.tracers.POP,
            NO₃ = model3.tracers.NO₃,
            Fe = model3.tracers.Fe)

simulation3.output_writers[:simple_output] =
            JLD2OutputWriter(model3, outputs, 
                            schedule = TimeInterval(1days), 
                            filename = "AMOC28_neweq",
                            overwrite_existing = true)

simulation3.output_writers[:checkpointer] = Checkpointer(model3,
            schedule = TimeInterval(new_eq_time),
            prefix = "AMOC28_checkpoint",
            overwrite_existing = false)

run!(simulation3, pickup = true)
=#
#################################### Video of all tracers ####################################

filepath = simulation2.output_writers[:simple_output].filepath
# filepath = "./AMOC25.jld2"

PO4_timeseries = FieldTimeSeries(filepath, "PO₄")
times = PO4_timeseries.times
xw, yw, zw = nodes(PO4_timeseries)
POP_timeseries = FieldTimeSeries(filepath, "POP")
DOP_timeseries = FieldTimeSeries(filepath, "DOP")
Fe_timeseries = FieldTimeSeries(filepath, "Fe")
NO3_timeseries = FieldTimeSeries(filepath, "NO₃")

n = Observable(1)
# title = @lift @sprintf("t = Year %d0", times[$n] / 3652.5days) 
title = @lift @sprintf("t = Day %d", times[$n] / 1days) 

# convert unit from mol/m³ to μM: 1e3*interior(...)
PO4ₙ = @lift 1e3*interior(PO4_timeseries[$n], 1, :, :)
avg_PO4ₙ = @lift mean(1e3*interior(PO4_timeseries[$n], 1, :, :), dims=1) 

POPₙ = @lift 1e3*interior(POP_timeseries[$n], 1, :, :)
avg_POPₙ = @lift mean(1e3*interior(POP_timeseries[$n], 1, :, :), dims=1) 

NO3ₙ = @lift 1e3*interior(NO3_timeseries[$n], 1, :, :)
Feₙ = @lift 1e6*interior(Fe_timeseries[$n], 1, :, :)
#=
fig = Figure(size = (1500, 1500))

ax_s = Axis(fig[2, 1]; xlabel = "y (km)", ylabel = "z (m)",title = "Stream function (m² s⁻¹)", aspect = 1)
hm_s = heatmap!(ax_s, yw./1e3, zw, Ψ; colormap = :viridis) 
Colorbar(fig[2, 2], hm_s; flipaxis = false)
contour!(ax_s, yw./1e3, grid.zᵃᵃᶜ[0:100], Ψ, levels = 10, color = :black)

ax_v = Axis(fig[2, 3]; xlabel = "y (km)", ylabel = "z (m)",title = "v (t) (cm s⁻¹)", aspect = 1)
hm_v = heatmap!(ax_v, yw./1e3, zw, vₙ; colorrange = (-3,3), colormap = :balance) 
Colorbar(fig[2, 4], hm_v; flipaxis = false)

ax_w = Axis(fig[2, 5]; xlabel = "y (km)", ylabel = "z (m)", title = "w (t) (cm s⁻¹)", aspect = 1)
hm_w = heatmap!(ax_w, yw./1e3, zw, wₙ; colorrange = (-8e-4, 8e-4), colormap = :balance) 
Colorbar(fig[2, 6], hm_w; flipaxis = false)

ax_PO4 = Axis(fig[3, 1]; xlabel = "y (km)", ylabel = "z (m)", title = "[PO₄] (μM)", aspect = 1)
hm_PO4 = heatmap!(ax_PO4, yw/1e3, zw, PO4ₙ; colorrange = (0,2.5),colormap = :rainbow1) 
Colorbar(fig[3, 2], hm_PO4; flipaxis = false)
contour!(ax_PO4, yw/1e3, zw[1:100], PO4ₙ, levels = 5, color = :black)

ax_NO3 = Axis(fig[3, 3]; xlabel = "y (km)", ylabel = "z (m)", title = "[NO₃] (μM)", aspect = 1)
hm_NO3 = heatmap!(ax_NO3, yw/1e3, zw, NO3ₙ; colorrange = (0,35),colormap = :rainbow1) 
Colorbar(fig[3, 4], hm_NO3; flipaxis = false)
contour!(ax_NO3, yw/1e3, zw[1:100], NO3ₙ, levels = 5, color = :black)

ax_Fe = Axis(fig[3, 5]; xlabel = "y (km)", ylabel = "z (m)", title = "[Fe] (nM)", aspect = 1)
hm_Fe = heatmap!(ax_Fe, yw/1e3, zw, Feₙ; colorrange = (0,1),colormap = :rainbow1) 
Colorbar(fig[3, 6], hm_Fe; flipaxis = false)
contour!(ax_Fe, yw/1e3, zw[1:100], Feₙ, levels = 5, color = :black)

ax_POP = Axis(fig[4, 3]; xlabel = "y (km)", ylabel = "z (m)", title = "[POP] (μM)", aspect = 1)
hm_POP = heatmap!(ax_POP, yw/1e3, zw, POPₙ; colorrange = (0,0.015),colormap = :rainbow1) 
Colorbar(fig[4, 4], hm_POP; flipaxis = false)

ax_avg_PO4 = Axis(fig[4, 1:2]; xlabel = "[PO₄] (μM)", ylabel = "z (m)", title = "Average [PO₄] (μM)", yaxisposition = :right)
xlims!(ax_avg_PO4, 0, 2)
PO4_prof = lines!(ax_avg_PO4, avg_PO4ₙ[][1, :], zw[1:100])

ax_avg_POP = Axis(fig[4, 5:6]; xlabel = "[POP] (μM)", ylabel = "z (m)", title = "Average [POP] (μM)",yaxisposition = :right)
xlims!(ax_avg_POP, 0, 0.015)
POP_prof = lines!(ax_avg_POP, avg_POPₙ[][1, :], zw[1:100])

fig[1, 1:6] = Label(fig, title, tellwidth=false)

# And, finally, we record a movie.
frames = 1:length(times)
record(fig, "AMOC28_perturb.mp4", frames, framerate=40) do i
    n[] = i
    PO4_prof[1] = avg_PO4ₙ[][1, :]
    POP_prof[1] = avg_POPₙ[][1, :]
end
nothing #hide
=#
############################## Compare initial vs final ##############################

# PO4_last = 1e3*interior(PO4_timeseries[end], 1, :, :) 
# POP_last = 1e3*interior(POP_timeseries[end], 1, :, :) 
# DOP_last = 1e3*interior(DOP_timeseries[end], 1, :, :) 
# sum(PO4_last.+POP_last.+DOP_last)

PO4_init = 1e3*interior(PO4_timeseries[1], 1, :, :) 
POP_init = 1e3*interior(POP_timeseries[1], 1, :, :) 
# DOP_init = 1e3*interior(DOP_timeseries[1], 1, :, :)
# sum(PO4_init.+POP_init.+DOP_init)

n = Observable(1)
title = @lift @sprintf("t = Day %d", times[$n] / 1days) 

diff_PO4 = @lift 1e3*(interior(PO4_timeseries[$n], 1, :, :) .- interior(PO4_timeseries[1], 1, :, :))
diff_POP = @lift 1e3*(interior(POP_timeseries[$n], 1, :, :) .- interior(POP_timeseries[1], 1, :, :))
# diff_DOP = @lift 1e3*(interior(DOP_timeseries[$n], 1, :, :) .- interior(DOP_timeseries[1], 1, :, :))
# sum(diff_PO4)

# One upwelling station
up_PO4ₙ = @lift 1e3*interior(PO4_timeseries[$n], 1, 10, :)
# up_DOPₙ = @lift 1e3*interior(DOP_timeseries[$n], 1, 10, :)
up_POPₙ = @lift 1e3*interior(POP_timeseries[$n], 1, 10, :)

init_up_PO4ₙ = 1e3*interior(PO4_timeseries[1], 1, 10, :)
init_up_POPₙ = 1e3*interior(POP_timeseries[1], 1, 10, :)
# init_up_DOPₙ = 1e3*interior(DOP_timeseries[1], 1, 10, :)

# One downwelling station
down_PO4ₙ = @lift 1e3*interior(PO4_timeseries[$n], 1, 490, :)
# down_DOPₙ = @lift 1e3*interior(DOP_timeseries[$n], 1, 490, :)
down_POPₙ = @lift 1e3*interior(POP_timeseries[$n], 1, 490, :)

init_down_PO4ₙ = 1e3*interior(PO4_timeseries[1], 1, 490, :)
init_down_POPₙ = 1e3*interior(POP_timeseries[1], 1, 490, :)
# init_down_DOPₙ = 1e3*interior(DOP_timeseries[1], 1, 490, :)

fig_compare = Figure(size=(900, 1100))

# Row 1: final - initial
ax_PO4 = Axis(fig_compare[2, 1]; xlabel = "y (km)", ylabel = "z (m)", title = "Δ[PO₄] (μM)", aspect = 1)
hm_PO4 = heatmap!(ax_PO4, yw/1e3, zw, diff_PO4; colorrange = (-0.15,0.15),colormap = :balance) 
Colorbar(fig_compare[2, 2], hm_PO4; flipaxis = false)

ax_POP = Axis(fig_compare[2, 3]; xlabel = "y (km)", ylabel = "z (m)", title = "Δ[POP] (μM)", aspect = 1)
hm_POP = heatmap!(ax_POP, yw/1e3, zw, diff_POP; colorrange = (-0.012,0.012),colormap = :balance) 
Colorbar(fig_compare[2, 4], hm_POP; flipaxis = false)

# ax_DOP = Axis(fig_compare[2, 5]; xlabel = "x (m)", ylabel = "z (m)", title = "Δ[DOP] (μM)", aspect = 1)
# hm_DOP = heatmap!(ax_DOP, yw, zw, diff_DOP; colorrange = (-0.06,0.06),colormap = :balance) 
# Colorbar(fig_compare[2, 6], hm_DOP; flipaxis = false)

# Row 2: upwelling site
up_PO4 = Axis(fig_compare[3, 1:2]; xlabel = "[PO₄] (μM)", ylabel = "z (m)", yaxisposition = :left)
xlims!(up_PO4, 2, 4)
lines!(up_PO4, init_up_PO4ₙ, zw[1:100],color = :red, linestyle = :dash, label = "Initial")
PO4_up = lines!(up_PO4, up_PO4ₙ[], zw[1:100], label = "Dynamic")
axislegend(up_PO4, position = :lb)

up_POP = Axis(fig_compare[3, 3:4]; xlabel = "[POP] (μM)", ylabel = "z (m)", yaxisposition = :left)
xlims!(up_POP, 0, 0.02)
lines!(up_POP, init_up_POPₙ, zw[1:100],color = :red, linestyle = :dash, label = "Initial")
POP_up = lines!(up_POP, up_POPₙ[], zw[1:100], label = "Dynamic")
axislegend(up_POP, position = :rb)

# up_DOP = Axis(fig_compare[3, 5:6]; xlabel = "[DOP] (μM)", ylabel = "z (m)", yaxisposition = :left)
# xlims!(up_DOP, -0.02, 0.5)
# lines!(up_DOP, init_up_DOPₙ, zw[1:200],color = :red, linestyle = :dash, label = "Initial")
# DOP_up = lines!(up_DOP, up_DOPₙ[], zw[1:200], label = "Dynamic")
# axislegend(up_DOP, position = :rb)

# Row 3: downwelling site
down_PO4 = Axis(fig_compare[4, 1:2]; xlabel = "[PO₄] (μM)", ylabel = "z (m)", yaxisposition = :left)
xlims!(down_PO4, 0, 1)
lines!(down_PO4, init_down_PO4ₙ, zw[1:100],color = :red, linestyle = :dash, label = "Initial")
PO4_down = lines!(down_PO4, down_PO4ₙ[], zw[1:100], label = "Dynamic")
axislegend(down_PO4, position = :lb)

down_POP = Axis(fig_compare[4, 3:4]; xlabel = "[POP] (μM)", ylabel = "z (m)", yaxisposition = :left)
xlims!(down_POP, 0, 0.003)
lines!(down_POP, init_down_POPₙ, zw[1:100],color = :red, linestyle = :dash, label = "Initial")
POP_down = lines!(down_POP, down_POPₙ[], zw[1:100], label = "Dynamic")
axislegend(down_POP, position = :rb)

# down_DOP = Axis(fig_compare[4, 5:6]; xlabel = "[DOP] (μM)", ylabel = "z (m)", yaxisposition = :left)
# xlims!(down_DOP, -0.01, 0.1)
# lines!(down_DOP, init_down_DOPₙ, zw[1:200],color = :red, linestyle = :dash, label = "Initial")
# DOP_down = lines!(down_DOP, down_DOPₙ[][1, :], zw[1:200], label = "Dynamic")
# axislegend(down_DOP, position = :rb)

fig_compare[1, 1:4] = Label(fig_compare, title, tellwidth=false)

# And, finally, we record a movie.
frames = 1:length(times)
record(fig_compare, "AMOC64_PARseasonal.mp4", frames, framerate=40) do i
    n[] = i
    PO4_up[1] = up_PO4ₙ[]#[1, :]
    POP_up[1] = up_POPₙ[]
    # DOP_up[1] = up_DOPₙ[]
    PO4_down[1] = down_PO4ₙ[]
    POP_down[1] = down_POPₙ[]
    # DOP_down[1] = down_DOPₙ[]
end
nothing #hide

################################# Plot fluxes ################################# 
#=
time_after = 365
 
POP_t_up = zeros(Nz,time_after) # zw,n
POP_t_up_diff = zeros(Nz,time_after) 

# Initial peak flux: POP[97] ~80 m
# Final peak flux: POP[94] ~140 m
F100 = zeros(1,time_after) 
F140 = zeros(1,time_after) 
F200 = zeros(1,time_after) 
F500 = zeros(1,time_after) 
F1000 = zeros(1,time_after) 
F2000 = zeros(1,time_after) 

for i in 1:time_after
    POP_t_up[:, i] = 1e3*interior(POP_timeseries[i], 1, 10, :)
    POP_t_up_diff[:, i] = 1e3*(interior(POP_timeseries[i], 1, 10, :).-interior(POP_timeseries[1], 1, 10, :))
    
    F100[1,i] = 10 .* POP_t_up[96,i] 
    F140[1,i] = 10 .* POP_t_up[94,i] 
    F200[1,i] = 10 .* POP_t_up[91,i] #/ POP_t_up[93,i] 
    F500[1,i] = 10 .* POP_t_up[76,i] 
    F1000[1,i] = 10 .* POP_t_up[51,i] 
    F2000[1,i] = 10 .* POP_t_up[1,i] 
end
z₀ = -log(0.01)*25 
F200_martin = ones(1,time_after) .* ((200+z₀)/(100+z₀))^-0.84
F500_martin = ones(1,time_after) .* ((500+z₀)/(100+z₀))^-0.84
F1000_martin = ones(1,time_after) .* ((1000+z₀)/(100+z₀))^-0.84
F2000_martin = ones(1,time_after) .* ((2000+z₀)/(100+z₀))^-0.84

fig_flux = Figure(size=(800, 600))

ax_F1000 = Axis(fig_flux[1, 1]; xlabel = "time (days)", ylabel = "F(z) (mmol m⁻² d⁻¹)", yaxisposition = :left)
lines!(ax_F1000, 1:time_after, F140[1,:], color = :orange, label = "140 m (final peak)")
lines!(ax_F1000, 1:time_after, F200[1,:], color = :black, label = "200 m")
lines!(ax_F1000, 1:time_after, F500[1,:], color = :green, label = "500 m")
lines!(ax_F1000, 1:time_after, F1000[1,:], color = :red, label = "1000 m")
lines!(ax_F1000, 1:time_after, F2000[1,:], color = :blue, label = "2000 m")
axislegend(ax_F1000, position = :rt)

ax_F140 = Axis(fig_flux[2, 1]; xlabel = "time (days)", ylabel = "F(z)/F₁₄₀", yaxisposition = :left)
ylims!(ax_F140, 0, 1)
lines!(ax_F140, 1:time_after, F200[1,:]./F140[1,:], color = :black, label = "200 m/140 m")
lines!(ax_F140, 1:time_after, F200_martin[1,:], color = :black, linestyle = :dash, label = "200 m/140 m (Martin)")
lines!(ax_F140, 1:time_after, F500[1,:]./F140[1,:], color = :green, label = "500 m/140 m")
lines!(ax_F140, 1:time_after, F500_martin[1,:], color = :green, linestyle = :dash, label = "500 m/140 m (Martin)")
lines!(ax_F140, 1:time_after, F1000[1,:]./F140[1,:], color = :red, label = "1000 m/140 m")
lines!(ax_F140, 1:time_after, F1000_martin[1,:], color = :red, linestyle = :dash, label = "1000 m/140 m (Martin)")
lines!(ax_F140, 1:time_after, F2000[1,:]./F140[1,:], color = :blue, label = "2000 m/140 m")
lines!(ax_F140, 1:time_after, F2000_martin[1,:], color = :blue, linestyle = :dash, label = "2000 m/140 m (Martin)")
axislegend(ax_F140, position = :rt)

display(fig_flux)
=#
