using GLMakie
using Printf
using Statistics

using Oceananigans
using Oceananigans.Units

filepath1 = "./P5_51y.jld2" 
filepath2 = "./P7_51y.jld2"
filepath3 = "./P6_51y.jld2"

NCP_timeseries1 = FieldTimeSeries(filepath1, "NCP")
POP_timeseries1 = FieldTimeSeries(filepath1, "POP")
Premin_timeseries1 = FieldTimeSeries(filepath1, "Premin")
times = Premin_timeseries1.times
xw, yw, zw = nodes(Premin_timeseries1)

NCP_timeseries2 = FieldTimeSeries(filepath2, "NCP")
Premin_timeseries2 = FieldTimeSeries(filepath2, "Premin")
POP_timeseries2 = FieldTimeSeries(filepath2, "POP")

NCP_timeseries3 = FieldTimeSeries(filepath3, "NCP")
Premin_timeseries3 = FieldTimeSeries(filepath3, "Premin")
POP_timeseries3 = FieldTimeSeries(filepath3, "POP")

#################################################################
##### Calculated remineralization: derived from Martin curve #####
#################################################################
 
dz = 20 # z grid gap
z₀ = log(0.01)*25 
Martin_factor = ((zw[190] .+ z₀) ./ (zw .+ z₀)).^0.84

######################################################################
########################## A single station ##########################
######################################################################
#=
POP_flux1 = zeros(365, 200) 
True_eratio1 = zeros(365, 1) 
True_Teff1 = zeros(365, 1)

for i in 1:365
    POP_flux1[i, :] = 1e4*interior(POP_timeseries1[i], 1, 50, :)
    True_eratio1[i, 1] = POP_flux1[i, 195]/sum(dz*1days*1e3*interior(NCP_timeseries1[i], 1, 50, :))
    True_Teff1[i, 1] = POP_flux1[i, 150]/sum(POP_flux1[i, 190])
end
Martin_value = ((zw[150]+z₀)/(zw[190]+z₀))^-0.84
Martin_Teff = repeat([Martin_value], 365, 1)
=#
#=
POP_ref1 = POP_flux1[:, 190]
POP_Martin1 = zeros(365,200) 
Martin_eratio1 = zeros(365, 1) 
Martin_Teff1 = zeros(365, 1) 

for i in 1:365
    temp_Martin1 = repeat([POP_ref1[i]], 1, 200)
    for k in 1:200
        temp_Martin1[:,k] .= temp_Martin1[:,k] .* Martin_factor[k]
    end
    POP_Martin1[i, :] = temp_Martin1
    # Calculate e-ratio (100 m/NCP)
    Martin_eratio1[i, 1] = POP_Martin1[i, 195]/sum(dz*1days*1e3*interior(NCP_timeseries1[i], 1, 50, :))
    # Calculate transfer efficiency (200 m/1000 m)
    Martin_Teff1[i, 1] = POP_Martin1[i, 150]/POP_Martin1[i, 190]
end
=#
######################################################################
########################## Domain average ##########################
######################################################################
#
POP_flux1 = zeros(365, 500, 200) 
True_eratio1 = zeros(365, 1) 
True_Teff1 = zeros(365, 1)

POP_flux2 = zeros(365, 500, 200) 
True_eratio2 = zeros(365, 1) 
True_Teff2 = zeros(365, 1) 

POP_flux3 = zeros(365, 500, 200) 
True_eratio3 = zeros(365, 1) 
True_Teff3 = zeros(365, 1) 

for i in 1:365
    # Compile 2D POP fluxes on each day
    POP_flux1[i, :, :] = 1e4*interior(POP_timeseries1[i], 1, :, :)
    POP_flux2[i, :, :] = 1e4*interior(POP_timeseries2[i], 1, :, :)
    POP_flux3[i, :, :] = 1e4*interior(POP_timeseries3[i], 1, :, :)
    # Calculate e-ratio (115 m/NCP)
    # True_eratio1[i, 1] = sum(POP_flux1[i, :, 195])/sum(dz*1days*1e3*interior(NCP_timeseries1[i], 1, :, :))
    # True_eratio2[i, 1] = sum(POP_flux2[i, :, 195])/sum(dz*1days*1e3*interior(NCP_timeseries2[i], 1, :, :))
    # True_eratio3[i, 1] = sum(POP_flux3[i, :, 190])/sum(dz*1days*1e3*interior(NCP_timeseries3[i], 1, :, :))
    # Calculate transfer efficiency (200 m/1000 m)
    True_Teff1[i, 1] = sum(POP_flux1[i, :, 150])/sum(POP_flux1[i, :, 190])
    True_Teff2[i, 1] = sum(POP_flux2[i, :, 150])/sum(POP_flux2[i, :, 190])
    True_Teff3[i, 1] = sum(POP_flux3[i, :, 150])/sum(POP_flux3[i, :, 190])
end
Martin_value = ((zw[150]+z₀)/(zw[190]+z₀))^-0.84
Martin_Teff = repeat([Martin_value], 365, 1)

#=
F1000_1 = sum(POP_flux2[:, :, 150];dims=1)
F200_1 = sum(POP_flux2[:, :, 190];dims=1)
mean_Teff1 = mean(F1000_1 ./ F200_1)

F115 = (sum(POP_flux2[:, :, 195];dims=1))[:]
F_NCP=zeros(365, 500)
for i in 1:365
    F_NCP[i,:]= sum(dz*1days*1e3*interior(NCP_timeseries2[i], 1, :, :);dims=2)
end
FP = sum(F_NCP;dims=1)[:]
meaneR = mean(F115 ./ FP)
=#

#=
POP_ref1 = POP_flux1[:, :, 190]
POP_ref2 = POP_flux2[:, :, 190]
POP_ref3 = POP_flux3[:, :, 190]
# POP flux in Martin curve (mmol m-2 d-1)
POP_Martin1 = zeros(365, 500, 200) 
POP_Martin2 = zeros(365, 500, 200) 
POP_Martin3 = zeros(365, 500, 200) 
Martin_eratio1 = zeros(365, 1) 
Martin_eratio2 = zeros(365, 1) 
Martin_eratio3 = zeros(365, 1) 
Martin_Teff1 = zeros(365, 1) 
Martin_Teff2 = zeros(365, 1) 
Martin_Teff3 = zeros(365, 1) 
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
    # Calculate e-ratio (100 m/NCP)
    Martin_eratio1[i, 1] = sum(POP_Martin1[i, :, 195])/sum(dz*1days*1e3*interior(NCP_timeseries1[i], 1, :, :))
    Martin_eratio2[i, 1] = sum(POP_Martin2[i, :, 195])/sum(dz*1days*1e3*interior(NCP_timeseries2[i], 1, :, :))
    Martin_eratio3[i, 1] = sum(POP_Martin3[i, :, 195])/sum(dz*1days*1e3*interior(NCP_timeseries3[i], 1, :, :))
    # Calculate transfer efficiency (200 m/1000 m)
    Martin_Teff1[i, 1] = sum(POP_Martin1[i, :, 150])/sum(POP_Martin1[i, :, 190])
    Martin_Teff2[i, 1] = sum(POP_Martin2[i, :, 150])/sum(POP_Martin2[i, :, 190])
    Martin_Teff3[i, 1] = sum(POP_Martin3[i, :, 150])/sum(POP_Martin3[i, :, 190])
end
=#

# t = 1:1:365
# MLD1 = 100 .+ 50 .* sinpi.(2 .* t ./ 365)
# MLD2 = 100 .+ 50 .* sinpi.(2 .* t ./ (365/4))
# MLD3 = 100 .+ 50 .* sinpi.(2 .* t ./ (365/12))

#################################################################
######################### Plot fluxes #########################
#################################################################

#=
fig = Figure(size=(1500, 1200))

############# e-ratio #############
ax_t_eR = Axis(fig[1, 1]; ylabel = "ratio", xlabel = "t (day)", title = "e-ratio (F₁₀₀/NCP)")
lines!(ax_t_eR, 1:1:365, vec(True_eratio1), linewidth = 2, label = "Model")
# lines!(ax_t_eR, 1:1:365,vec(Martin_eratio1),linewidth = 2, label = "Martin")
axislegend(ax_t_eR, position = :rb)
ylims!(ax_t_eR, 0.15, 0.25)

ax_t_MLD = Axis(fig[1, 1]; ylabel = "(m)", xlabel = "t (day)", yaxisposition=:right, 
                yticks=0:50:200, yticklabelcolor=:red3, ylabelcolor=:red3)
lines!(ax_t_MLD, 1:1:365, MLD1, linewidth = 1, color=:red3, label = "MLD")
ylims!(ax_t_MLD, 200, 0)
linkxaxes!(ax_t_eR, ax_t_MLD)
axislegend(ax_t_MLD, position = :rt)

ax_t_eR2 = Axis(fig[1, 2]; ylabel = "ratio", xlabel = "t (day)", title = "e-ratio (F₁₀₀/NCP)")
lines!(ax_t_eR2, 1:1:365, vec(True_eratio2), linewidth = 2, label = "Model")
# lines!(ax_t_eR2, 1:1:365,vec(Martin_eratio2),linewidth = 2, label = "Martin")
axislegend(ax_t_eR2, position = :rb)
ylims!(ax_t_eR2, 0.15, 0.25)
ax_t_MLD2 = Axis(fig[1, 2]; ylabel = "(m)", xlabel = "t (day)", yaxisposition=:right, 
                yticks=0:50:200, yticklabelcolor=:red3, ylabelcolor=:red3)
lines!(ax_t_MLD2, 1:1:365, MLD2, linewidth = 1, color=:red3, label = "MLD")
ylims!(ax_t_MLD2, 200,0)
linkxaxes!(ax_t_eR2, ax_t_MLD2)
axislegend(ax_t_MLD2, position = :rt)

ax_t_eR3 = Axis(fig[1, 3]; ylabel = "ratio", xlabel = "t (day)", title = "e-ratio (F₁₀₀/NCP)")
lines!(ax_t_eR3, 1:1:365, vec(True_eratio3), linewidth = 2, label = "Model")
# lines!(ax_t_eR3, 1:1:365,vec(Martin_eratio3),linewidth = 2, label = "Martin")
axislegend(ax_t_eR3, position = :rb)
ylims!(ax_t_eR3, 0.15, 0.25)
ax_t_MLD3 = Axis(fig[1, 3]; ylabel = "(m)", xlabel = "t (day)", yaxisposition=:right, 
                yticks=0:50:200, yticklabelcolor=:red3, ylabelcolor=:red3)
lines!(ax_t_MLD3, 1:1:365, MLD3, linewidth = 1, color=:red3, label = "MLD")
ylims!(ax_t_MLD3, 200, 0)
linkxaxes!(ax_t_eR3, ax_t_MLD3)
axislegend(ax_t_MLD3, position = :rt)

############# Teff #############
ax_t_tot = Axis(fig[2, 1]; ylabel = "ratio", xlabel = "t (day)", title = "Transfer effeciency (F₁₀₀₀/F₂₀₀)")
lines!(ax_t_tot, 1:1:365, vec(True_Teff1), linewidth = 2, label = "Model")
lines!(ax_t_tot, 1:1:365,vec(Martin_Teff1),linewidth = 2, label = "Martin")
axislegend(ax_t_tot, position = :rb)
ylims!(ax_t_tot, 0.3, 0.4)
ax_MLD = Axis(fig[2, 1]; ylabel = "(m)", xlabel = "t (day)", yaxisposition=:right, 
                yticks=0:50:200, yticklabelcolor=:red3, ylabelcolor=:red3)
lines!(ax_MLD, 1:1:365, MLD1, linewidth = 1, color=:red3, label = "MLD")
ylims!(ax_MLD, 200, 0)
linkxaxes!(ax_t_tot, ax_MLD)
axislegend(ax_MLD, position = :lt)

ax_t_tot2 = Axis(fig[2, 2]; ylabel = "ratio", xlabel = "t (day)", title = "Transfer effeciency (F₁₀₀₀/F₂₀₀)")
lines!(ax_t_tot2, 1:1:365, vec(True_Teff2), linewidth = 2, label = "Model")
lines!(ax_t_tot2, 1:1:365,vec(Martin_Teff2),linewidth = 2, label = "Martin")
axislegend(ax_t_tot2, position = :rb)
ylims!(ax_t_tot2, 0.3, 0.4)
ax_MLD2 = Axis(fig[2, 2]; ylabel = "(m)", xlabel = "t (day)", yaxisposition=:right, 
                yticks=0:50:200, yticklabelcolor=:red3, ylabelcolor=:red3)
lines!(ax_MLD2, 1:1:365, MLD2, linewidth = 1, color=:red3, label = "MLD")
ylims!(ax_MLD2, 200, 0)
linkxaxes!(ax_t_tot2, ax_MLD2)
axislegend(ax_MLD2, position = :rt)

ax_t_tot3 = Axis(fig[2, 3]; ylabel = "ratio", xlabel = "t (day)", title = "Transfer effeciency (F₁₀₀₀/F₂₀₀)")
lines!(ax_t_tot3, 1:1:365, vec(True_Teff3), linewidth = 2, label = "Model")
lines!(ax_t_tot3, 1:1:365,vec(Martin_Teff3),linewidth = 2, label = "Martin")
axislegend(ax_t_tot3, position = :rb)
ylims!(ax_t_tot3, 0.3, 0.4)
ax_MLD3 = Axis(fig[2, 3]; ylabel = "(m)", xlabel = "t (day)", yaxisposition=:right, 
                yticks=0:50:200, yticklabelcolor=:red3, ylabelcolor=:red3)
lines!(ax_MLD3, 1:1:365, MLD3, linewidth = 1, color=:red3, label = "MLD")
ylims!(ax_MLD3, 200, 0)
linkxaxes!(ax_t_tot3, ax_MLD3)
axislegend(ax_MLD3, position = :rt)

display(fig)

=#
# function deviation_from_average(vec)
#     avg = mean(vec)
#     max_dev = (maximum(vec) - avg) / avg * 100
#     min_dev = (minimum(vec) - avg) / avg * 100
#     return max_dev, min_dev
# end
# max_dev, min_dev = deviation_from_average(vec(True_eratio1))
#################################################################
######################### Plot Teff together #########################
#################################################################
#
fig = Figure(size=(600, 500))

# ax_t_er = Axis(fig[1, 1]; ylabel = "ratio", xlabel = "t (day)", title = "e-ratio (F₁₁₅/NCP)")
# lines!(ax_t_er, 1:1:365, vec(True_eratio1), linewidth = 2, label = "Year")
# lines!(ax_t_er, 1:1:365, vec(True_eratio2), linewidth = 2, label = "Year+Month")
# lines!(ax_t_er, 1:1:365, vec(True_eratio3), linewidth = 2, label = "Year+Week")
# ylims!(ax_t_er, 0,0.5)

ax_t_tot = Axis(fig[1, 1]; ylabel = "ratio", xlabel = "t (day)", title = "Transfer effeciency (F₁₀₀₀/F₂₀₀)")
lines!(ax_t_tot, 1:1:365, vec(True_Teff1), linewidth = 2, label = "T=365")
lines!(ax_t_tot, 1:1:365, vec(True_Teff2), linewidth = 2, label = "T=160")
lines!(ax_t_tot, 1:1:365, vec(True_Teff3), linewidth = 2, label = "T=30")
lines!(ax_t_tot, 1:1:365,vec(Martin_Teff),linewidth = 2, label = "Martin")
axislegend(ax_t_tot, position = :lt)
ylims!(ax_t_tot, 0, 1.3)

display(fig)
# Mavg = mean(vec(Martin_Teff1)) #0.35246585432250294
#
