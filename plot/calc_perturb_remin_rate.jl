using GLMakie
using Printf
using Statistics

using Oceananigans
using Oceananigans.Units
using Oceananigans.Fields: ZeroField, CenterField, FunctionField
using Oceananigans.BoundaryConditions: fill_halo_regions!

#################################################################
##### Load data: Premin (rate) and POP concentration (flux) #####
#################################################################

filepath1 = "./AMOC115_1perY.jld2"
filepath2 = "./AMOC115_4perY.jld2"
filepath3 = "./AMOC115_12perY.jld2"

Premin_timeseries1 = FieldTimeSeries(filepath1, "Premin")
times = Premin_timeseries1.times
xw, yw, zw = nodes(Premin_timeseries1)
Premin_timeseries2 = FieldTimeSeries(filepath2, "Premin")
Premin_timeseries3 = FieldTimeSeries(filepath3, "Premin")

POP_timeseries1 = FieldTimeSeries(filepath1, "POP")
POP_timeseries2 = FieldTimeSeries(filepath2, "POP")
POP_timeseries3 = FieldTimeSeries(filepath3, "POP")

#################################################################
######### True remineralization: direct model outputs #########
#################################################################

avg_Premin1 = zeros(365, 200) # (time, z_grids)
avg_Premin2 = zeros(365, 200) 
avg_Premin3 = zeros(365, 200) 
tot_Premin1 = zeros(365, 1)
tot_Premin2 = zeros(365, 1)
tot_Premin3 = zeros(365, 1)

for i in 1:365
    avg_Premin1[i,:] = mean(1days*1e3*interior(Premin_timeseries1[i], 1, :, :), dims=1) 
    avg_Premin2[i,:] = mean(1days*1e3*interior(Premin_timeseries2[i], 1, :, :), dims=1)
    avg_Premin3[i,:] = mean(1days*1e3*interior(Premin_timeseries3[i], 1, :, :), dims=1)
    tot_Premin1[i,1] = sum(1days*1e3*interior(Premin_timeseries1[i], 1, :, 150:190)) 
    tot_Premin2[i,1] = sum(1days*1e3*interior(Premin_timeseries2[i], 1, :, 150:190))  
    tot_Premin3[i,1] = sum(1days*1e3*interior(Premin_timeseries3[i], 1, :, 150:190))  
end

#################################################################
##### Calculated remineralization: derived from Martin curve #####
#################################################################
 
z₀ = log(0.01)*25 
Martin_factor = ((zw[190] .+ z₀) ./ (zw .+ z₀)).^0.84

POP_flux1 = zeros(365, 500, 200) 
POP_flux2 = zeros(365, 500, 200) 
POP_flux3 = zeros(365, 500, 200) 
for i in 1:365
    POP_flux1[i, :, :] = 1e4*interior(POP_timeseries1[i], 1, :, :)
    POP_flux2[i, :, :] = 1e4*interior(POP_timeseries2[i], 1, :, :)
    POP_flux3[i, :, :] = 1e4*interior(POP_timeseries3[i], 1, :, :)
end

POP_ref1 = POP_flux1[:, :, 190]
POP_ref2 = POP_flux2[:, :, 190]
POP_ref3 = POP_flux3[:, :, 190]
# POP flux in Martin curve (mmol m-2 d-1)
POP_Martin1 = zeros(365, 500, 200) 
POP_Martin2 = zeros(365, 500, 200) 
POP_Martin3 = zeros(365, 500, 200) 
for i in 1:365
    temp_Martin1 = repeat(POP_ref1[i,:], 1, 200)
    temp_Martin2 = repeat(POP_ref2[i,:], 1, 200)
    temp_Martin3 = repeat(POP_ref3[i,:], 1, 200)
    for k in 1:200
        temp_Martin1[:,k] .= temp_Martin1[:,k] .* Martin_factor[k]
        temp_Martin2[:,k] .= temp_Martin2[:,k] .* Martin_factor[k]
        temp_Martin3[:,k] .= temp_Martin3[:,k] .* Martin_factor[k]
    end
    POP_Martin1[i, :, :] = temp_Martin1
    POP_Martin2[i, :, :] = temp_Martin2
    POP_Martin3[i, :, :] = temp_Martin3
end

# Based on the Martin POP flux, calculate remineralization rates
dz = -20.0  # depth intervals = 20 m
w = 10.0 # Sinking rate = 10 m/day

Martin_R1 = zeros(365,500,200)
Martin_R2 = zeros(365,500,200)
Martin_R3 = zeros(365,500,200)
for t in 1:365
    # Compute central differences
    for i in 1:500
        for k in 2:199
            Martin_R1[t, i, k] = -(POP_Martin1[t, i, k+1] - POP_Martin1[t, i, k-1]) / (2 * dz)
            Martin_R2[t, i, k] = -(POP_Martin2[t, i, k+1] - POP_Martin2[t, i, k-1]) / (2 * dz)
            Martin_R3[t, i, k] = -(POP_Martin3[t, i, k+1] - POP_Martin3[t, i, k-1]) / (2 * dz)
        end
    end
    # Forward and backward differences for boundaries
    for i in 1:500
        Martin_R1[t, i, 1] = -(POP_Martin1[t, i, 2] - POP_Martin1[t, i, 1]) / dz
        Martin_R1[t, i, 200] = -(POP_Martin1[t, i, 200] - POP_Martin1[t, i, 199]) / dz
        Martin_R2[t, i, 1] = -(POP_Martin2[t, i, 2] - POP_Martin2[t, i, 1]) / dz
        Martin_R2[t, i, 200] = -(POP_Martin2[t, i, 200] - POP_Martin2[t, i, 199]) / dz
        Martin_R3[t, i, 1] = -(POP_Martin3[t, i, 2] - POP_Martin3[t, i, 1]) / dz
        Martin_R3[t, i, 200] = -(POP_Martin3[t, i, 200] - POP_Martin3[t, i, 199]) / dz
    end
end

avg_Martin1 = zeros(365, 200) # (time, z_grids)
avg_Martin2 = zeros(365, 200) 
avg_Martin3 = zeros(365, 200) 
tot_Martin1 = zeros(365, 1)
tot_Martin2 = zeros(365, 1)
tot_Martin3 = zeros(365, 1)
for t in 1:365
    avg_Martin1[t,:] = mean(Martin_R1[t,:,:], dims=1)
    avg_Martin2[t,:] = mean(Martin_R2[t,:,:], dims=1)
    avg_Martin3[t,:] = mean(Martin_R3[t,:,:], dims=1)
    tot_Martin1[t,1] = sum(Martin_R1[t,:,150:190])  
    tot_Martin2[t,1] = sum(Martin_R2[t,:,150:190]) 
    tot_Martin3[t,1] = sum(Martin_R3[t,:,150:190]) 
end

year_Martin1 = sum(tot_Martin1) # mmol P m⁻³ y⁻¹
year_Martin2 = sum(tot_Martin2) 
year_Martin3 = sum(tot_Martin3) 

#################################################################
######################### Plot! #########################
#################################################################

fig = Figure(size=(1500, 1200))

############# total #############
ax_t_tot = Axis(fig[1, 1:2]; ylabel = "POP remin rate (mmol m⁻³ d⁻¹)", xlabel = "t (day)", title = "Total remin at 200 - 1000 m")
lines!(ax_t_tot, 1:1:365, vec(tot_Premin1), linewidth = 2, label = "True")
lines!(ax_t_tot, 1:1:365,vec(tot_Martin1),linewidth = 2, label = "Martin")
axislegend(ax_t_tot, position = :rt)
ylims!(ax_t_tot, 0.7, 0.91)

ax_t_tot2 = Axis(fig[1, 3:4]; ylabel = "POP remin rate (mmol m⁻³ d⁻¹)", xlabel = "t (day)", title = "Total remin at 200 - 1000 m")
lines!(ax_t_tot2, 1:1:365, vec(tot_Premin2), linewidth = 2, label = "True")
lines!(ax_t_tot2, 1:1:365, vec(tot_Martin2),linewidth = 2, label = "Martin")
axislegend(ax_t_tot2, position = :rb)
ylims!(ax_t_tot2, 0.7, 0.91)

ax_t_tot3 = Axis(fig[1, 5:6]; ylabel = "POP remin rate (mmol m⁻³ d⁻¹)", xlabel = "t (day)", title = "Total remin at 200 - 1000 m")
lines!(ax_t_tot3, 1:1:365, vec(tot_Premin3), linewidth = 2, label = "True")
lines!(ax_t_tot3, 1:1:365, vec(tot_Martin3),linewidth = 2, label = "Martin")
axislegend(ax_t_tot3, position = :rb)
ylims!(ax_t_tot3, 0.7, 0.91)

############# True heatmap #############
ax_True1 = Axis(fig[2, 1]; xlabel = "t (days)", ylabel = "z (m)", title = "Ave. POP remin rate: log₁₀(mmol m⁻³ d⁻¹)", aspect = 1)
hm_True1 = heatmap!(ax_True1, 1:1:365, zw, log10.(avg_Premin1); colorrange = (-4.9, -3.9),colormap = :rainbow1) 
Colorbar(fig[2, 2], hm_True1; flipaxis = false)
ylims!(ax_True1, -1000, -200)

ax_True2 = Axis(fig[2, 3]; xlabel = "t (days)", ylabel = "z (m)", title = "Ave. POP remin rate: log₁₀(mmol m⁻³ d⁻¹)", aspect = 1)
hm_True2 = heatmap!(ax_True2, 1:1:365, zw, log10.(avg_Premin2); colorrange = (-4.9, -3.9),colormap = :rainbow1) 
Colorbar(fig[2, 4], hm_True2; flipaxis = false)
ylims!(ax_True2, -1000, -200)

ax_True3 = Axis(fig[2, 5]; xlabel = "t (days)", ylabel = "z (m)", title = "Ave. POP remin rate: log₁₀(mmol m⁻³ d⁻¹)", aspect = 1)
hm_True3 = heatmap!(ax_True3, 1:1:365, zw, log10.(avg_Premin3); colorrange = (-4.9, -3.9),colormap = :rainbow1) 
Colorbar(fig[2, 6], hm_True3; flipaxis = false)
ylims!(ax_True3, -1000, -200)

############# Martin heatmap #############
ax_Martin1 = Axis(fig[3, 1]; xlabel = "t (days)", ylabel = "z (m)", title = "Ave. Martin-derived remin rate: log₁₀(mmol m⁻³ d⁻¹)", aspect = 1)
hm_Martin1 = heatmap!(ax_Martin1, 1:1:365, zw, log10.(avg_Martin1); colorrange = (-4.9, -3.9),colormap = :rainbow1) 
Colorbar(fig[3, 2], hm_Martin1; flipaxis = false)
ylims!(ax_Martin1, -1000, -200)

ax_Martin2 = Axis(fig[3, 3]; xlabel = "t (days)", ylabel = "z (m)", title = "Ave. Martin-derived remin rate: log₁₀(mmol m⁻³ d⁻¹)", aspect = 1)
hm_Martin2 = heatmap!(ax_Martin2, 1:1:365, zw, log10.(avg_Martin2); colorrange = (-4.9, -3.9),colormap = :rainbow1) 
Colorbar(fig[3, 4], hm_Martin2; flipaxis = false)
ylims!(ax_Martin2, -1000, -200)

ax_Martin3 = Axis(fig[3, 5]; xlabel = "t (days)", ylabel = "z (m)", title = "Ave. Martin-derived remin rate: log₁₀(mmol m⁻³ d⁻¹)", aspect = 1)
hm_Martin3 = heatmap!(ax_Martin3, 1:1:365, zw, log10.(avg_Martin3); colorrange = (-4.9, -3.9),colormap = :rainbow1) 
Colorbar(fig[3, 6], hm_Martin3; flipaxis = false)
ylims!(ax_Martin3, -1000, -200)

############# Diff heatmap #############
ax_diff1 = Axis(fig[4, 1]; xlabel = "t (days)", ylabel = "z (m)", title = "True - Martin-derived remin rate (mmol m⁻³ d⁻¹)", aspect = 1)
hm_diff1 = heatmap!(ax_diff1, 1:1:365, zw, (avg_Premin1 .- avg_Martin1); colorrange = (-6e-6,6e-6),colormap = :balance) 
Colorbar(fig[4, 2], hm_diff1; flipaxis = false)
ylims!(ax_diff1, -1000, -200)

ax_diff2 = Axis(fig[4, 3]; xlabel = "t (days)", ylabel = "z (m)", title = "True - Martin-derived remin rate (mmol m⁻³ d⁻¹)", aspect = 1)
hm_diff2 = heatmap!(ax_diff2, 1:1:365, zw, (avg_Premin2 .- avg_Martin2); colorrange = (-6e-6,6e-6),colormap = :balance) 
Colorbar(fig[4, 4], hm_diff2; flipaxis = false)
ylims!(ax_diff2, -1000, -200)

ax_diff3 = Axis(fig[4, 5]; xlabel = "t (days)", ylabel = "z (m)", title = "True - Martin-derived remin rate (mmol m⁻³ d⁻¹)", aspect = 1)
hm_diff3 = heatmap!(ax_diff3, 1:1:365, zw, (avg_Premin3 .- avg_Martin3); colorrange = (-6e-6,6e-6),colormap = :balance) 
Colorbar(fig[4, 6], hm_diff3; flipaxis = false)
ylims!(ax_diff3, -1000, -200)

display(fig)