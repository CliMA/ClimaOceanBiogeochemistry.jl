# Plot vertical profiles of POP fluxes: true vs Martin curve
# Focus: different time slices -> deviation?

using GLMakie
using Printf
using Statistics

using Oceananigans
using Oceananigans.Units
using XLSX, JLD2, DataFrames

#################################################################
##### Load data: Premin (rate) and POP concentration (flux) #####
#################################################################

filepath1 = "./P5_51y.jld2"
filepath2 = "./P6_51y.jld2"

POP_timeseries1 = FieldTimeSeries(filepath1, "POP")
POP_timeseries2 = FieldTimeSeries(filepath2, "POP")
times = POP_timeseries1.times
xw, yw, zw = nodes(POP_timeseries1)

#################################################################
##### Calculated remineralization: derived from Martin curve #####
#################################################################
 
z₀ = log(0.01)*25 
Martin_factor = ((zw[191] + z₀) ./ (zw .+ z₀)).^0.84
PtoCratio = 117

POC_flux1 = zeros(365, 200) 
POC_flux2 = zeros(365, 200) 
for i in 1:365
    POC_flux1[i, :, :] = PtoCratio*1e4*interior(POP_timeseries1[i], 1, 115, :)
    POC_flux2[i, :, :] = PtoCratio*1e4*interior(POP_timeseries2[i], 1, 115, :)
end

POC_ref1 = POC_flux1[:, 191]
POC_ref2 = POC_flux2[:, 191]

# POC flux in Martin curve (mmol m-2 d-1)
POC_Martin1 = zeros(365, 200) 
POC_Martin2 = zeros(365, 200) 

const_r = 0.017
exp_factor = exp.(const_r/10 .* (zw.-z₀))
POC_Exponential1 = zeros(365, 200) 
POC_Exponential2 = zeros(365, 200) 
ref1 = POC_flux1[:, 195]
ref2 = POC_flux2[:, 195]

for i in 1:365
    POC_Martin1[i, :] = repeat([POC_ref1[i]], 1, 200).* Martin_factor'
    POC_Martin2[i, :] = repeat([POC_ref2[i]], 1, 200).* Martin_factor'

    POC_Exponential1[i, :] = repeat([ref1[i]], 1, 200).* exp_factor'
    POC_Exponential2[i, :] = repeat([ref2[i]], 1, 200).* exp_factor'
end

#################################################################
###### Pick out time slices of POC flux: True vs. Martin ######
#################################################################

# Annual sum
# Ftrue1 = vec(sum(POC_flux1; dims=1))
# Fmartin1 = vec(sum(POC_Martin1; dims=1)) 

# Ftrue2 = vec(sum(POC_flux2; dims=1))
# Fmartin2 = vec(sum(POC_Martin2; dims=1)) 

# Annual mean
Ftrue1 = vec(mean(POC_flux1; dims=1))
Ftrue1_std = vec(std(POC_flux1; dims=1))
Fmartin1 = vec(mean(POC_Martin1; dims=1)) 
Fmartin1_std = vec(std(POC_Martin1; dims=1)) 
Fex1 = vec(mean(POC_Exponential1; dims=1)) 
Fex1_std = vec(std(POC_Exponential1; dims=1)) 

Ftrue2 = vec(mean(POC_flux2; dims=1))
Ftrue2_std = vec(std(POC_flux2; dims=1))
Fmartin2 = vec(mean(POC_Martin2; dims=1)) 
Fmartin2_std = vec(std(POC_Martin2; dims=1)) 
Fex2 = vec(mean(POC_Exponential2; dims=1)) 
Fex2_std = vec(std(POC_Exponential2; dims=1)) 

# day 100
Ftrue1_100 = POC_flux1[95,:]
Fmartin1_100 = POC_Martin1[95,:]
Fex1_100 = POC_Exponential1[95,:]

Ftrue2_100 = POC_flux2[95,:]
Fmartin2_100 = POC_Martin2[95,:]
Fex2_100 = POC_Exponential2[95,:]

# day 180
Ftrue1_180 = POC_flux1[180,:]
Fmartin1_180 = POC_Martin1[180,:]
Fex1_180 = POC_Exponential1[180,:]

Ftrue2_180 = POC_flux2[180,:]
Fmartin2_180 = POC_Martin2[180,:]
Fex2_180 = POC_Exponential2[180,:]

# day 315
Ftrue1_315 = POC_flux1[320,:]
Fmartin1_315 = POC_Martin1[320,:]
Fex1_315 = POC_Exponential1[320,:]

Ftrue2_315 = POC_flux2[320,:]
Fmartin2_315 = POC_Martin2[320,:]
Fex2_315 = POC_Exponential2[320,:]

#################################################################
######################### Load obs data #########################
#################################################################
#
@load "data/POCflux_Cael2018.jld2" obs_POCflux
obs_depth_Cael = obs_POCflux[:,1]
obs_POC_flux_Cael = obs_POCflux[:,2]./ 12  # mg C m⁻² d⁻¹ to mmol C m⁻² d⁻¹
@load "data/POCflux_Cael2018Martin.jld2" obs_POCflux
obs_depth_Martin = obs_POCflux[:,1]
obs_POC_flux_Martin = obs_POCflux[:,2]./ 12 # mg C m⁻² d⁻¹ to mol C m⁻² y⁻¹

#################################################################
######################### Plot fluxes #########################
#################################################################
#
zw_band = vcat(zw, reverse(zw))
Ftrue1_band = vcat(Ftrue1 .+ Ftrue1_std, reverse(Ftrue1 .- Ftrue1_std))
Ftrue2_band = vcat(Ftrue2 .+ Ftrue2_std, reverse(Ftrue2 .- Ftrue2_std))

# Fmartin1_band = vcat(Fmartin1 .+ Fmartin1_std, reverse(Fmartin1 .- Fmartin1_std))
# zw_band_adj = vcat(zw_band[150:190], zw_band[210:250])
# Fmartin1_band_adj = vcat(Fmartin1_band[150:190],Fmartin1_band[210:250])

#
fig = Figure(size=(500, 500))
### Annual total ###
ax1 = Axis(fig[1, 1]; ylabel = "Depth (m)", xlabel = "POC flux (mmol m⁻² d⁻¹)", title = "POC flux (daily average)")
scatter!(ax1, obs_POC_flux_Cael, -obs_depth_Cael, marker = :circle, color = (:grey66,0.5), markersize=6, label = "Observations")
lines!(ax1, Ftrue1, zw, color=:dodgerblue4, linewidth = 10, label = "Yearly ± 1 std")
poly!(Ftrue1_band, zw_band, color = (:dodgerblue4, 0.2)) 
lines!(ax1, Ftrue2, zw, color=:plum3, linewidth = 5, label = "Monthly ± 1 std")
poly!(Ftrue2_band, zw_band, color = (:plum3, 0.3)) 
lines!(ax1, Fmartin1[150:192], zw[150:192], color=:seagreen, linewidth = 3, label = "Martin curve")
lines!(ax1, Fex1[150:192], zw[150:192], color=:sienna4, linewidth = 3, linestyle=:dash, label = "Exponential")
ylims!(ax1, -1000, 0)
xlims!(ax1, 0, 15)
axislegend(ax1, position = :rb)

display(fig)
#
#=
fig = Figure(size=(1000, 500))
############# Day 100 #############
ax11 = Axis(fig[1, 1]; ylabel = "Depth (m)", xlabel = "POC flux (mmol m⁻² d⁻¹)", title = "POC flux (Day 95)")
lines!(ax11, Ftrue1_100, zw, color=:dodgerblue4, linewidth = 3, label = "Model")
lines!(ax11, Fmartin1_100[150:195], zw[150:195], color=:seagreen, linewidth = 2, label = "Martin curve")
lines!(ax11, Fex1_100[150:195], zw[150:195], color=:sienna4, linewidth = 2, linestyle=:dash, label = "Exponential")
ylims!(ax11, -1000, 0)
xlims!(ax11, 0, 15)
axislegend(ax11, position = :rb)

ax21 = Axis(fig[2, 1]; ylabel = "Depth (m)", xlabel = "POC flux (mmol m⁻² d⁻¹)", title = "POC flux (Day 95)")
lines!(ax21, Ftrue2_100, zw, color=:plum3, linewidth = 3, label = "Model")
lines!(ax21, Fmartin2_100[150:195], zw[150:195], color=:seagreen, linewidth = 2, label = "Martin")
lines!(ax21, Fex2_100[150:195], zw[150:195], color=:sienna4, linewidth = 2, linestyle=:dash, label = "Exponential")
ylims!(ax21, -1000, 0)
xlims!(ax21, 0, 15)
axislegend(ax21, position = :rb)

############# Day 180 #############
ax12 = Axis(fig[1, 2]; ylabel = "Depth (m)", xlabel = "POC flux (mmol m⁻² d⁻¹)", title = "POC flux (Day 180)")
lines!(ax12, Ftrue1_180, zw, color=:dodgerblue4, linewidth = 3, label = "Model")
lines!(ax12, Fmartin1_180[150:195], zw[150:195], color=:seagreen, linewidth = 2, label = "Martin")
lines!(ax12, Fex1_180[150:195], zw[150:195], color=:sienna4, linewidth = 2,linestyle=:dash,  label = "Exponential")
ylims!(ax12, -1000, 0)
xlims!(ax12, 0, 15)

ax22 = Axis(fig[2, 2]; ylabel = "Depth (m)", xlabel = "POC flux (mmol m⁻² d⁻¹)", title = "POC flux (Day 180)")
lines!(ax22, Ftrue2_180, zw, color=:plum3, linewidth = 3, label = "Model")
lines!(ax22, Fmartin2_180[150:195], zw[150:195], color=:seagreen, linewidth = 2, label = "Martin")
lines!(ax22, Fex2_180[150:195], zw[150:195], color=:sienna4, linewidth = 2, linestyle=:dash, label = "Exponential")
ylims!(ax22, -1000, 0)
xlims!(ax22, 0, 15)

############# Day 315 #############
ax13 = Axis(fig[1, 3]; ylabel = "Depth (m)", xlabel = "POC flux (mmol m⁻² d⁻¹)", title = "POC flux (Day 320)")
lines!(ax13, Ftrue1_315, zw, color=:dodgerblue4, linewidth = 3, label = "Model")
lines!(ax13, Fmartin1_315[150:195], zw[150:195], color=:seagreen, linewidth = 2, label = "Martin")
lines!(ax13, Fex1_315[150:195], zw[150:195], color=:sienna4, linewidth = 2, linestyle=:dash, label = "Exponential")
ylims!(ax13, -1000, 0)
xlims!(ax13, 0, 15)

ax23 = Axis(fig[2, 3]; ylabel = "Depth (m)", xlabel = "POC flux (mmol m⁻² d⁻¹)", title = "POC flux (Day 320)")
lines!(ax23, Ftrue2_315, zw, color=:plum3, linewidth = 3, label = "Model")
lines!(ax23, Fmartin2_315[150:195], zw[150:195], color=:seagreen, linewidth = 2, label = "Martin")
lines!(ax23, Fex2_315[150:195], zw[150:195], color=:sienna4, linewidth = 2, linestyle=:dash, label = "Exponential")
ylims!(ax23, -1000, 0)
xlims!(ax23, 0, 15)

display(fig)
=#