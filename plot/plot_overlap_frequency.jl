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

filepath1 = "./P4_21y.jld2"

NCP_timeseries1 = FieldTimeSeries(filepath1, "NCP")
Premin_timeseries1 = FieldTimeSeries(filepath1, "Premin")
times = Premin_timeseries1.times
xw, yw, zw = nodes(Premin_timeseries1)
Dremin_timeseries1 = FieldTimeSeries(filepath1, "Dremin")

POP_timeseries1 = FieldTimeSeries(filepath1, "POP")

#################################################################
######### True remineralization: direct model outputs #########
#################################################################

tot_NCP1 = zeros(365, 1)
tot_remin1 = zeros(365, 1)

avg_Premin1 = zeros(365, 200) # (time, z_grids)
tot_Premin1 = zeros(365, 1)

for i in 1:365
    tot_NCP1[i,1] = sum(1days*1e3*interior(NCP_timeseries1[i], 1, :, :)) 
    tot_remin1[i,1] = sum(1days*1e3*interior(Premin_timeseries1[i], 1, :, :)) .+ sum(1days*1e3*interior(Dremin_timeseries1[i], 1, :, :))

    avg_Premin1[i,:] = mean(1days*1e3*interior(Premin_timeseries1[i], 1, :, :), dims=1) 
    tot_Premin1[i,1] = sum(1days*1e3*interior(Premin_timeseries1[i], 1, :, 150:190)) 
end

#################################################################
##### Calculated remineralization: derived from Martin curve #####
#################################################################
 
z₀ = log(0.01)*25 
Martin_factor = ((zw[190] .+ z₀) ./ (zw .+ z₀)).^0.84

POP_flux1 = zeros(365, 500, 200) 

avg_TrueFlux1 = zeros(365, 200) 
True_Teff1 = zeros(365, 1) 
for i in 1:365
    # Compile 2D POP fluxes on each day
    POP_flux1[i, :, :] = 1e4*interior(POP_timeseries1[i], 1, :, :)
    # Calculate domain-average flux 
    avg_TrueFlux1[i,:] = mean(1e4*interior(POP_timeseries1[i], 1, :, :), dims=1) 
    # Calculate transfer efficiency (200 m/1000 m)
    True_Teff1[i, 1] = sum(POP_flux1[i, :, 150])/sum(POP_flux1[i, :, 190])
end

POP_ref1 = POP_flux1[:, :, 190]
# POP flux in Martin curve (mmol m-2 d-1)
POP_Martin1 = zeros(365, 500, 200)  
avg_MartinFlux1 = zeros(365, 200) 
Martin_Teff1 = zeros(365, 1) 
for i in 1:365
    temp_Martin1 = repeat(POP_ref1[i,:], 1, 200)
    for k in 1:200
        temp_Martin1[:,k] .= temp_Martin1[:,k] .* Martin_factor[k]
    end
    POP_Martin1[i, :, :] = temp_Martin1
    # Calculate domain-average flux 
    avg_MartinFlux1[i,:] = mean(POP_Martin1[i, :, :], dims=1) 
    # Calculate transfer efficiency (200 m/1000 m)
    Martin_Teff1[i, 1] = sum(POP_Martin1[i, :, 150])/sum(POP_Martin1[i, :, 190])
end

# Based on the Martin POP flux, calculate remineralization rates 
# rP = -(1/w)*(dP/dz) = -dF/dz
dz = -20.0  # depth intervals = 20 m

Martin_R1 = zeros(365,500,200)
for t in 1:365
    # Compute central differences
    for i in 1:500
        for k in 2:199
            Martin_R1[t, i, k] = -(POP_Martin1[t, i, k+1] - POP_Martin1[t, i, k-1]) / (2 * dz)
        end
    end
    # Forward and backward differences for boundaries
    for i in 1:500
        Martin_R1[t, i, 1] = -(POP_Martin1[t, i, 2] - POP_Martin1[t, i, 1]) / dz
        Martin_R1[t, i, 200] = -(POP_Martin1[t, i, 200] - POP_Martin1[t, i, 199]) / dz
    end
end

avg_Martin1 = zeros(365, 200) # (time, z_grids)
tot_Martin1 = zeros(365, 1)
for t in 1:365
    avg_Martin1[t,:] = mean(Martin_R1[t,:,:], dims=1)
    tot_Martin1[t,1] = sum(Martin_R1[t,:,150:190])  
end

year_Martin1 = sum(tot_Martin1) # mmol P m⁻³ y⁻¹

#################################################################
######################### Plot fluxes #########################
#################################################################

fig = Figure(size=(1100, 800))

############# Transfer efficiency #############
ax_t_Teff = Axis(fig[1, 1:2]; ylabel = "ratio", xlabel = "t (day)", title = "Transfer effeciency (F₁₀₀₀/F₂₀₀)")
lines!(ax_t_Teff, 1:1:365, vec(True_Teff1), linewidth = 2, label = "True")
lines!(ax_t_Teff, 1:1:365,vec(Martin_Teff1),linewidth = 2, label = "Martin")
axislegend(ax_t_Teff, position = :rb)
ylims!(ax_t_Teff, 0.25, 0.45)

############# Total remin #############
ax_t_totR = Axis(fig[1, 3:4]; ylabel = "POP remin rate (mmol m⁻³ d⁻¹)", xlabel = "t (day)", title = "Total remin at 200 - 1000 m")
lines!(ax_t_totR, 1:1:365, vec(tot_Premin1), linewidth = 2, label = "True")
lines!(ax_t_totR, 1:1:365,vec(tot_Martin1),linewidth = 2, label = "Martin")
axislegend(ax_t_totR, position = :rt)
ylims!(ax_t_totR, 0.7, 1.1)

ax_t_totP = Axis(fig[1, 5:6]; ylabel = "(mmol m⁻³ d⁻¹)", xlabel = "t (day)", title = "Integrated NCP & Remin")
lines!(ax_t_totP, 1:1:365, vec(tot_NCP1), linewidth = 2, label = "Production")
lines!(ax_t_totP, 1:1:365,vec(tot_remin1),linewidth = 2, label = "Remineralization")
axislegend(ax_t_totP, position = :lb)
ylims!(ax_t_totP, 6, 9)

############# True heatmap #############
ax_True1 = Axis(fig[2, 1]; xlabel = "t (days)", ylabel = "z (m)", title = "Ave. POP flux (mmol m⁻² d⁻¹)", aspect = 1)
hm_True1 = heatmap!(ax_True1, 1:1:365, zw, avg_TrueFlux1; colorrange = (0.01, 0.06),colormap = :thermal) 
Colorbar(fig[2, 2], hm_True1; flipaxis = false)
ylims!(ax_True1, -1000, -200)

############# Martin heatmap #############
ax_Martin1 = Axis(fig[2, 3]; xlabel = "t (days)", ylabel = "z (m)", title = "Ave. Martin-derived flux (mmol m⁻² d⁻¹)", aspect = 1)
hm_Martin1 = heatmap!(ax_Martin1, 1:1:365, zw, avg_MartinFlux1; colorrange = (0.01, 0.06),colormap = :thermal) 
Colorbar(fig[2, 4], hm_Martin1; flipaxis = false)
ylims!(ax_Martin1, -1000, -200)

############# Diff heatmap #############
ax_diff1 = Axis(fig[2, 5]; xlabel = "t (days)", ylabel = "z (m)", title = "True - Martin-derived flux (mmol m⁻² d⁻¹)", aspect = 1)
hm_diff1 = heatmap!(ax_diff1, 1:1:365, zw, (avg_TrueFlux1 .- avg_MartinFlux1); colorrange = (-0.006, 0.006),colormap = :balance) 
Colorbar(fig[2, 6], hm_diff1; flipaxis = false)
ylims!(ax_diff1, -1000, -200)

display(fig)