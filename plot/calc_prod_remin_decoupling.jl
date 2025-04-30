using GLMakie
using Printf
using Statistics

using Oceananigans
using Oceananigans.Units

filepath1 = "./P5_51y.jld2"
filepath2 = "./P6_51y.jld2"
# filepath3 = "./P3_21y.jld2"

#################################################################
##### Load data: rates of NCP, Dremin, Premin #####
#################################################################

NCP_timeseries1 = FieldTimeSeries(filepath1, "NCP")
Dremin_timeseries1 = FieldTimeSeries(filepath1, "Dremin")
Premin_timeseries1 = FieldTimeSeries(filepath1, "Premin")

NCP_timeseries2 = FieldTimeSeries(filepath2, "NCP")
Dremin_timeseries2 = FieldTimeSeries(filepath2, "Dremin")
Premin_timeseries2 = FieldTimeSeries(filepath2, "Premin")

# NCP_timeseries3 = FieldTimeSeries(filepath3, "NCP")
# Dremin_timeseries3 = FieldTimeSeries(filepath3, "Dremin")
# Premin_timeseries3 = FieldTimeSeries(filepath3, "Premin")

times = Premin_timeseries1.times
xw, yw, zw = nodes(Premin_timeseries1)

#################################################################
######### prod vs. remin #########
#################################################################

avg_NCP1 = zeros(365, 200) # (time, z_grids)
avg_NCP2 = zeros(365, 200) 
# avg_NCP3 = zeros(365, 200) 
tot_NCP1 = zeros(365, 1)
tot_NCP2 = zeros(365, 1)
# tot_NCP3 = zeros(365, 1)

avg_Dremin1 = zeros(365, 200) # (time, z_grids)
avg_Dremin2 = zeros(365, 200) 
# avg_Dremin3 = zeros(365, 200) 
tot_Dremin1 = zeros(365, 1)
tot_Dremin2 = zeros(365, 1)
# tot_Dremin3 = zeros(365, 1)

avg_Premin1 = zeros(365, 200) # (time, z_grids)
avg_Premin2 = zeros(365, 200) 
# avg_Premin3 = zeros(365, 200) 
tot_Premin1 = zeros(365, 1)
tot_Premin2 = zeros(365, 1)
# tot_Premin3 = zeros(365, 1)

for i in 1:365
    # NCP
    avg_NCP1[i,:] = mean(1days*1e3*interior(NCP_timeseries1[i], 1, :, :), dims=1) 
    avg_NCP2[i,:] = mean(1days*1e3*interior(NCP_timeseries2[i], 1, :, :), dims=1)
    # avg_NCP3[i,:] = mean(1days*1e3*interior(NCP_timeseries3[i], 1, :, :), dims=1)
    tot_NCP1[i,1] = sum(1days*1e3*interior(NCP_timeseries1[i], 1, :, :)) 
    tot_NCP2[i,1] = sum(1days*1e3*interior(NCP_timeseries2[i], 1, :, :))  
    # tot_NCP3[i,1] = sum(1days*1e3*interior(NCP_timeseries3[i], 1, :, :))  
    # DOP remin
    avg_Dremin1[i,:] = mean(1days*1e3*interior(Dremin_timeseries1[i], 1, :, :), dims=1) 
    avg_Dremin2[i,:] = mean(1days*1e3*interior(Dremin_timeseries2[i], 1, :, :), dims=1)
    # avg_Dremin3[i,:] = mean(1days*1e3*interior(Dremin_timeseries3[i], 1, :, :), dims=1)
    tot_Dremin1[i,1] = sum(1days*1e3*interior(Dremin_timeseries1[i], 1, :, :)) 
    tot_Dremin2[i,1] = sum(1days*1e3*interior(Dremin_timeseries2[i], 1, :, :))  
    # tot_Dremin3[i,1] = sum(1days*1e3*interior(Dremin_timeseries3[i], 1, :, :))  
    # POP remin
    avg_Premin1[i,:] = mean(1days*1e3*interior(Premin_timeseries1[i], 1, :, :), dims=1) 
    avg_Premin2[i,:] = mean(1days*1e3*interior(Premin_timeseries2[i], 1, :, :), dims=1)
    # avg_Premin3[i,:] = mean(1days*1e3*interior(Premin_timeseries3[i], 1, :, :), dims=1)
    tot_Premin1[i,1] = sum(1days*1e3*interior(Premin_timeseries1[i], 1, :, :)) 
    tot_Premin2[i,1] = sum(1days*1e3*interior(Premin_timeseries2[i], 1, :, :))  
    # tot_Premin3[i,1] = sum(1days*1e3*interior(Premin_timeseries3[i], 1, :, :))  
end

avg_remin1 = avg_Dremin1 .+ avg_Premin1
tot_remin1 = tot_Dremin1 .+ tot_Premin1
avg_PoverR1 = avg_NCP1 ./ avg_remin1 

avg_remin2 = avg_Dremin2 .+ avg_Premin2
avg_PoverR2 = avg_NCP2 ./ avg_remin2 
tot_remin2 = tot_Dremin2 .+ tot_Premin2

# avg_remin3 = avg_Dremin3 .+ avg_Premin3
# avg_PoverR3 = avg_NCP3 ./ avg_remin3 
# tot_remin3 = tot_Dremin3 .+ tot_Premin3

# t = 1:1:365
# MLD1 = 100 .+ 90 .* sinpi.(2 .* t ./ 365)
# MLD2 = 100 .+ 90 .* sinpi.(2 .* t ./ (365/30))
# MLD3 = 100 .+ 50 .* sinpi.(2 .* t ./ (365/12))

#################################################################
######################### Plot fluxes #########################
#################################################################

fig = Figure(size=(800, 300))
############# integrated with time #############
ax_t_tot = Axis(fig[1, 1]; ylabel = "(mmol m⁻³ d⁻¹)", xlabel = "t (day)", title = "Integrated NCP & Remin")
lines!(ax_t_tot, 1:1:365, vec(tot_NCP1), linewidth = 3, label = "Production")
lines!(ax_t_tot, 1:1:365,vec(tot_remin1),linewidth = 3, label = "Total remineralization")
lines!(ax_t_tot, 1:1:365,vec(tot_Premin1),linewidth = 3, label = "POP remineralization")

axislegend(ax_t_tot, position = :rt)
ylims!(ax_t_tot, 0, 12)
# argmax(tot_NCP1) # find the index of the peak value
# ax_t_MLD = Axis(fig[1, 1]; ylabel = "(m)", xlabel = "t (day)", yaxisposition=:right, 
#                 yticks=0:50:200, yticklabelcolor=:red3, ylabelcolor=:red3)
# lines!(ax_t_MLD, 1:1:365, MLD1, linewidth = 1, color=:red3, label = "MLD")
# ylims!(ax_t_MLD, 200,0)
# linkxaxes!(ax_t_tot, ax_t_MLD)

ax_t_tot2 = Axis(fig[1, 2]; ylabel = "(mmol m⁻³ d⁻¹)", xlabel = "t (day)", title = "Integrated NCP & Remin")
lines!(ax_t_tot2, 1:1:365, vec(tot_NCP2), linewidth = 3, label = "Production")
lines!(ax_t_tot2, 1:1:365,vec(tot_remin2),linewidth = 3, label = "Total remineralization")
lines!(ax_t_tot2, 1:1:365,vec(tot_Premin2),linewidth = 3, label = "POP remineralization")

# axislegend(ax_t_tot2, position = :lb)
ylims!(ax_t_tot2, 0,12)
# ax_t_MLD2 = Axis(fig[2,1]; ylabel = "(m)", xlabel = "t (day)", yaxisposition=:right, 
#                 yticks=0:50:200, yticklabelcolor=:red3, ylabelcolor=:red3)
# lines!(ax_t_MLD2, 1:1:365, MLD2, linewidth = 1, color=:red3, label = "MLD")
# ylims!(ax_t_MLD2, 200,0)
# linkxaxes!(ax_t_tot2, ax_t_MLD2)

# ax_t_tot3 = Axis(fig[1, 3]; ylabel = "(mmol m⁻³ d⁻¹)", xlabel = "t (day)", title = "Integrated NCP & Remin")
# lines!(ax_t_tot3, 1:1:365, vec(tot_NCP3), linewidth = 3, label = "Production")
# lines!(ax_t_tot3, 1:1:365,vec(tot_remin3),linewidth = 3, label = "Remineralization")
# axislegend(ax_t_tot3, position = :lb)
# ylims!(ax_t_tot3, 5.8, 8.2)
# ax_t_MLD3 = Axis(fig[1, 5:6]; ylabel = "(m)", xlabel = "t (day)", yaxisposition=:right, 
#                 yticks=0:50:200, yticklabelcolor=:red3, ylabelcolor=:red3)
# lines!(ax_t_MLD3, 1:1:365, MLD3, linewidth = 1, color=:red3, label = "MLD")
# ylims!(ax_t_MLD3, 200,0)
# linkxaxes!(ax_t_tot3, ax_t_MLD3)

display(fig)

############# heatmap #############
# ax_ratio1 = Axis(fig[2, 1]; xlabel = "t (days)", ylabel = "z (m)", title = "Ave. ratio of NCP/remin", aspect = 1)
# hm_ratio1 = heatmap!(ax_ratio1, 1:1:365, zw, avg_PoverR1; colorrange = (0, 2),colormap = :bam)
# Colorbar(fig[2, 2], hm_ratio1; flipaxis = false)
# ylims!(ax_ratio1, -200, 0)
# contour!(ax_ratio1, 1:1:365, zw, avg_PoverR1, levels = [1], color = :black, linewidth=2)
# lines!(ax_ratio1, 1:1:365, -MLD1, linewidth = 2,color = :red3)

# ax_ratio2 = Axis(fig[2, 3]; xlabel = "t (days)", ylabel = "z (m)", title = "Ave. ratio of NCP/remin", aspect = 1)
# hm_ratio2 = heatmap!(ax_ratio2, 1:1:365, zw, avg_PoverR2; colorrange = (0, 2),colormap = :bam) 
# Colorbar(fig[2, 4], hm_ratio2; flipaxis = false)
# ylims!(ax_ratio2, -200, 0)
# contour!(ax_ratio2, 1:1:365, zw, avg_PoverR2, levels = [1], color = :black, linewidth=2)
# lines!(ax_ratio2, 1:1:365, -MLD2, linewidth = 2,color = :red3)

# ax_ratio3 = Axis(fig[2, 5]; xlabel = "t (days)", ylabel = "z (m)", title = "Ave. ratio of NCP/remin", aspect = 1)
# hm_ratio3 = heatmap!(ax_ratio3, 1:1:365, zw, avg_PoverR3; colorrange = (0, 2),colormap = :bam) 
# Colorbar(fig[2, 6], hm_ratio3; flipaxis = false)
# ylims!(ax_ratio3, -200, 0)
# contour!(ax_ratio3, 1:1:365, zw, avg_PoverR3, levels = [1], color = :black, linewidth=2)
# lines!(ax_ratio3, 1:1:365, -MLD3, linewidth = 2,color = :red3)

# display(fig)
############# specific time points #############
#=
ax_t_r = Axis(fig[3, 1:2]; ylabel = "z (m)", xlabel = "ratio", title = "NCP : Remin")
lines!(ax_t_r, avg_PoverR1[1,:], zw, linewidth = 2, label = "t = Day 1")
lines!(ax_t_r, avg_PoverR1[91,:], zw, linewidth = 2, label = "t = Day 91")
lines!(ax_t_r, avg_PoverR1[182,:], zw, linewidth = 2, label = "t = Day 182")
lines!(ax_t_r, avg_PoverR1[274,:], zw, linewidth = 2, label = "t = Day 274")
lines!(ax_t_r, avg_PoverR1[365,:], zw, linewidth = 2, label = "t = Day 365")
axislegend(ax_t_r, position = :rb)
ylims!(ax_t_r, -200, 0)
xlims!(ax_t_r, 0, 3)

ax_t_r2 = Axis(fig[3, 3:4]; ylabel = "z (m)", xlabel = "ratio", title = "NCP : Remin")
lines!(ax_t_r2, avg_PoverR2[1,:], zw, linewidth = 2, label = "t = Day 1")
lines!(ax_t_r2, avg_PoverR2[91,:], zw, linewidth = 2, label = "t = Day 91")
lines!(ax_t_r2, avg_PoverR2[182,:], zw, linewidth = 2, label = "t = Day 182")
lines!(ax_t_r2, avg_PoverR2[274,:], zw, linewidth = 2, label = "t = Day 274")
lines!(ax_t_r2, avg_PoverR2[365,:], zw, linewidth = 2, label = "t = Day 365")
axislegend(ax_t_r2, position = :rb)
ylims!(ax_t_r2, -200, 0)
xlims!(ax_t_r2, 0, 3)

ax_t_r3 = Axis(fig[3, 5:6]; ylabel = "z (m)", xlabel = "ratio", title = "NCP : Remin")
lines!(ax_t_r3, avg_PoverR3[1,:], zw, linewidth = 2, label = "t = Day 1")
lines!(ax_t_r3, avg_PoverR3[91,:], zw, linewidth = 2, label = "t = Day 91")
lines!(ax_t_r3, avg_PoverR3[182,:], zw, linewidth = 2, label = "t = Day 182")
lines!(ax_t_r3, avg_PoverR3[274,:], zw, linewidth = 2, label = "t = Day 274")
lines!(ax_t_r3, avg_PoverR3[365,:], zw, linewidth = 2, label = "t = Day 365")
axislegend(ax_t_r3, position = :rb)
ylims!(ax_t_r3, -200, 0)
xlims!(ax_t_r3, 0, 3)
=#
