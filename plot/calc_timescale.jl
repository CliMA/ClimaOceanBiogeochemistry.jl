using GLMakie
using Printf
using Statistics

using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients
using Oceananigans
using Oceananigans.Units

filepath1 = "./AMOC115_auxiliary.jld2"

#################################################################
##### Load data: rates of NCP, Dremin, Premin #####
#################################################################

NCP_timeseries = FieldTimeSeries(filepath1, "NCP")
DOP_timeseries = FieldTimeSeries(filepath1, "DOP")
POP_timeseries = FieldTimeSeries(filepath1, "POP")

times = NCP_timeseries.times
xw, yw, zw = nodes(NCP_timeseries)

# mol P m-3 d-1
NCP_baseline = 1days*interior(NCP_timeseries[end], 1, :, 190:200) 
# mol P m-3
POP_baseline = interior(POP_timeseries[end], 1, :, 190:200) 
DOP_baseline = interior(DOP_timeseries[end], 1, :, 190:200)

# days
NCP_timescale = (POP_baseline .+ DOP_baseline)./NCP_baseline
avg_NCP_timescale = mean(NCP_timescale; dims = 1)

# POP remin 
z₀ = log(0.01)*25
r = -0.84*10 ./ (zw .+ z₀)
Premin_timescale = 1 ./ r

# DOP remin 
Dremin = 1/(2/365.25)
Dremin_timescale = Dremin .* ones(200)

#################################################################
######################### Plot timescales #########################
#################################################################

fig = Figure(size=(1000, 1000))

ax_NCP = Axis(fig[1, 1]; xlabel = "y (10³ km)", ylabel = "z (m)", title = "NCP timescale (days)", aspect = 1)
hm_NCP = heatmap!(ax_NCP, yw/1e6, zw[190:200], NCP_timescale; colorrange = (0, 200), colormap = :RdYlGn_6)  
Colorbar(fig[1, 2], hm_NCP; flipaxis = false)
ylims!(ax_NCP, -150, 0)
contour!(ax_NCP, yw/1e6, zw[190:200], NCP_timescale, levels = [91], color = :black, linewidth=2)

ax_avgNCP = Axis(fig[1, 3:4]; xlabel = "days", ylabel = "z (m)", title = "mean NCP timescale", aspect = 1)
lines!(ax_avgNCP, vec(avg_NCP_timescale),zw[190:200], linewidth = 2)
ylims!(ax_avgNCP, -150, 0)
xlims!(ax_avgNCP, 0, 500)

ax_Premin = Axis(fig[2, 1:2]; xlabel = "days", ylabel = "z (m)", title = "POP remin timescale", aspect = 1)
lines!(ax_Premin, Premin_timescale[150:200],zw[150:200], linewidth = 2)
ylims!(ax_Premin, -1000, 0)
xlims!(ax_Premin, 0, 150)

ax_Dremin = Axis(fig[2, 3:4]; xlabel = "days", ylabel = "z (m)", title = "DOP remin timescale", aspect = 1)
lines!(ax_Dremin, Dremin_timescale[150:200],zw[150:200], linewidth = 2)
xlims!(ax_Dremin, 170, 190)
ylims!(ax_Dremin, -1000, 0)

display(fig)