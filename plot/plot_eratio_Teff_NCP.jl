using GLMakie
using Printf
using Statistics

using Oceananigans
using Oceananigans.Units

filepath1 = "./P5_51y.jld2"
filepath2 = "./P6_51y.jld2"
filepath3 = "./P7_51y.jld2"

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
########################## Domain average ##########################
######################################################################
#
POP_flux1 = zeros(365, 500, 200) 
True_eratio1 = zeros(365, 1) 
True_Teff1 = zeros(365, 1)
# True_T100_1 = zeros(365, 1)

POP_flux2 = zeros(365, 500, 200) 
True_eratio2 = zeros(365, 1) 
True_Teff2 = zeros(365, 1) 
# True_T100_2 = zeros(365, 1)

tot_NCP1 = zeros(365, 1)
tot_NCP2 = zeros(365, 1)
tot_NCP3 = zeros(365, 1)

POP_flux3 = zeros(365, 500, 200) 
True_eratio3 = zeros(365, 1) 
True_Teff3 = zeros(365, 1) 

for i in 1:365
    # Compile 2D POP fluxes on each day
    POP_flux1[i, :, :] = 1e4*interior(POP_timeseries1[i], 1, :, :)
    POP_flux2[i, :, :] = 1e4*interior(POP_timeseries2[i], 1, :, :)
    POP_flux3[i, :, :] = 1e4*interior(POP_timeseries3[i], 1, :, :)
    # Calculate NCP
    tot_NCP1[i,1] = sum(1days*1e3*interior(NCP_timeseries1[i], 1, :, :)) 
    tot_NCP2[i,1] = sum(1days*1e3*interior(NCP_timeseries2[i], 1, :, :))  
    tot_NCP3[i,1] = sum(1days*1e3*interior(NCP_timeseries3[i], 1, :, :)) 
    # Calculate e-ratio (115 m/NCP)
    True_eratio1[i, 1] = sum(POP_flux1[i, :, 195])/sum(dz*1days*1e3*interior(NCP_timeseries1[i], 1, :, :))
    True_eratio2[i, 1] = sum(POP_flux2[i, :, 195])/sum(dz*1days*1e3*interior(NCP_timeseries2[i], 1, :, :))
    True_eratio3[i, 1] = sum(POP_flux3[i, :, 190])/sum(dz*1days*1e3*interior(NCP_timeseries3[i], 1, :, :))
    # Calculate T100 (215 m/115 m)
    # True_T100_1[i, 1] = sum(POP_flux1[i, :, 190])/sum(POP_flux1[i, :, 195])
    # True_T100_2[i, 1] = sum(POP_flux2[i, :, 190])/sum(POP_flux2[i, :, 195])
    # Calculate transfer efficiency (200 m/1000 m)
    True_Teff1[i, 1] = sum(POP_flux1[i, :, 150])/sum(POP_flux1[i, :, 190])
    True_Teff2[i, 1] = sum(POP_flux2[i, :, 150])/sum(POP_flux2[i, :, 190])
    True_Teff3[i, 1] = sum(POP_flux3[i, :, 150])/sum(POP_flux3[i, :, 190])
end
Martin_value = ((zw[150]+z₀)/(zw[190]+z₀))^-0.84
Martin_Teff = repeat([Martin_value], 365, 1)

#################################################################
######################### Plot fluxes #########################
#################################################################

#
fig = Figure(size=(1500, 1500))

############# NCP vs. e-ratio #############
ax_t_eR = Axis(fig[1, 1]; ylabel = "e-ratio (F₁₁₅/NCP)", xlabel = "t (day)", title = "yearly: e-ratio vs. NCP",
            yticklabelcolor=:brown3, ylabelcolor=:brown3)
lines!(ax_t_eR, 1:1:365, vec(True_eratio1), color=:brown3, linewidth = 2, label = "e-ratio")
axislegend(ax_t_eR, position = :lt)
ylims!(ax_t_eR, 0, 0.5)

ax_t_MLD = Axis(fig[1, 1]; ylabel = "NCP (mmol m⁻³ d⁻¹)", xlabel = "t (day)", yaxisposition=:right, 
               yticklabelcolor=:royalblue1, ylabelcolor=:royalblue1)
lines!(ax_t_MLD, 1:1:365, vec(tot_NCP1), linewidth = 2, color=:royalblue1, label = "NCP")
ylims!(ax_t_MLD,  2,12)
linkxaxes!(ax_t_eR, ax_t_MLD)
axislegend(ax_t_MLD, position = :rt)

ax_t_eR2 = Axis(fig[1, 2]; ylabel = "e-ratio (F₁₁₅/NCP)", xlabel = "t (day)", title = "max Δϕ: e-ratio vs. NCP",
            yticklabelcolor=:brown3, ylabelcolor=:brown3)
lines!(ax_t_eR2, 1:1:365, vec(True_eratio3), color=:brown3, linewidth = 2, label = "e-ratio")
axislegend(ax_t_eR2, position = :lt)
ylims!(ax_t_eR2, 0, 0.5)
ax_t_MLD2 = Axis(fig[1, 2]; ylabel = "NCP (mmol m⁻³ d⁻¹)", xlabel = "t (day)", yaxisposition=:right, 
                yticklabelcolor=:royalblue1, ylabelcolor=:royalblue1)
lines!(ax_t_MLD2, 1:1:365, vec(tot_NCP3), linewidth = 2, color=:royalblue1, label = "NCP")
ylims!(ax_t_MLD2, 2,12)
linkxaxes!(ax_t_eR2, ax_t_MLD2)
axislegend(ax_t_MLD2, position = :rt)

ax_t_eR3 = Axis(fig[1, 3]; ylabel = "e-ratio (F₁₁₅/NCP)", xlabel = "t (day)", title = "monthly: e-ratio vs. NCP",
            yticklabelcolor=:brown3, ylabelcolor=:brown3)
lines!(ax_t_eR3, 1:1:365, vec(True_eratio2), color=:brown3, linewidth = 2, label = "e-ratio")
axislegend(ax_t_eR3, position = :lt)
ylims!(ax_t_eR3, 0, 0.5)
ax_t_MLD3 = Axis(fig[1, 3]; ylabel = "NCP (mmol m⁻³ d⁻¹)", xlabel = "t (day)", yaxisposition=:right, 
                yticklabelcolor=:royalblue1, ylabelcolor=:royalblue1)
lines!(ax_t_MLD3, 1:1:365, vec(tot_NCP2), linewidth = 2, color=:royalblue1, label = "NCP")
ylims!(ax_t_MLD3, 2,12)
linkxaxes!(ax_t_eR3, ax_t_MLD3)
axislegend(ax_t_MLD3, position = :rt)
#
############# NCP vs. T100 #############
#=
ax_t_eR = Axis(fig[2, 1]; ylabel = "T₁₀₀ (F₂₁₅/F₁₁₅)", xlabel = "t (day)", title = "yearly",
            yticklabelcolor=:brown3, ylabelcolor=:brown3)
lines!(ax_t_eR, 1:1:365, vec(True_T100_1), color=:brown3, linewidth = 2, label = "T₁₀₀")
axislegend(ax_t_eR, position = :lt)
ylims!(ax_t_eR, 0.5, 1.5)

ax_t_MLD = Axis(fig[2, 1]; ylabel = "NCP (mmol m⁻³ d⁻¹)", xlabel = "t (day)", yaxisposition=:right, 
               yticklabelcolor=:royalblue1, ylabelcolor=:royalblue1)
lines!(ax_t_MLD, 1:1:365, vec(tot_NCP1), linewidth = 2, color=:royalblue1, label = "NCP")
ylims!(ax_t_MLD,  2,12)
linkxaxes!(ax_t_eR, ax_t_MLD)
axislegend(ax_t_MLD, position = :rt)

ax_t_eR2 = Axis(fig[2, 2]; ylabel = "T₁₀₀ (F₂₁₅/F₁₁₅)", xlabel = "t (day)", title = "monthly",
            yticklabelcolor=:brown3, ylabelcolor=:brown3)
lines!(ax_t_eR2, 1:1:365, vec(True_T100_2), color=:brown3, linewidth = 2, label = "T₁₀₀")
axislegend(ax_t_eR2, position = :lt)
ylims!(ax_t_eR2, 0.5, 1.5)
ax_t_MLD2 = Axis(fig[2, 2]; ylabel = "NCP (mmol m⁻³ d⁻¹)", xlabel = "t (day)", yaxisposition=:right, 
                yticklabelcolor=:royalblue1, ylabelcolor=:royalblue1)
lines!(ax_t_MLD2, 1:1:365, vec(tot_NCP2), linewidth = 2, color=:royalblue1, label = "NCP")
ylims!(ax_t_MLD2, 2,12)
linkxaxes!(ax_t_eR2, ax_t_MLD2)
axislegend(ax_t_MLD2, position = :rt)
=#
############# NCP vs. Teff #############
ax_t_eR = Axis(fig[2, 1]; ylabel = "Tₑ (F₁₀₀₀/F₂₀₀)", xlabel = "t (day)", title = "yearly: Tₑ vs. NCP",
            yticklabelcolor=:darkgoldenrod2, ylabelcolor=:darkgoldenrod2)
lines!(ax_t_eR, 1:1:365, vec(True_Teff1), color=:darkgoldenrod2, linewidth = 2, label = "Tₑ")
axislegend(ax_t_eR, position = :lt)
ylims!(ax_t_eR, 0, 1.3)
ax_t_MLD = Axis(fig[2, 1]; ylabel = "NCP (mmol m⁻³ d⁻¹)", xlabel = "t (day)", yaxisposition=:right, 
               yticklabelcolor=:royalblue1, ylabelcolor=:royalblue1)
lines!(ax_t_MLD, 1:1:365, vec(tot_NCP1), linewidth = 2, color=:royalblue1, label = "NCP")
ylims!(ax_t_MLD,  2,12)
linkxaxes!(ax_t_eR, ax_t_MLD)
axislegend(ax_t_MLD, position = :rt)

ax_t_eR2 = Axis(fig[2, 2]; ylabel = "Tₑ (F₁₀₀₀/F₂₀₀)", xlabel = "t (day)", title = "max Δϕ: Tₑ vs. NCP",
            yticklabelcolor=:darkgoldenrod2, ylabelcolor=:darkgoldenrod2)
lines!(ax_t_eR2, 1:1:365, vec(True_Teff3), color=:darkgoldenrod2, linewidth = 2, label = "Tₑ")
axislegend(ax_t_eR2, position = :lt)
ylims!(ax_t_eR2, 0, 1.3)
ax_t_MLD2 = Axis(fig[2, 2]; ylabel = "NCP (mmol m⁻³ d⁻¹)", xlabel = "t (day)", yaxisposition=:right, 
                yticklabelcolor=:royalblue1, ylabelcolor=:royalblue1)
lines!(ax_t_MLD2, 1:1:365, vec(tot_NCP3), linewidth = 2, color=:royalblue1, label = "NCP")
ylims!(ax_t_MLD2, 2,12)
linkxaxes!(ax_t_eR2, ax_t_MLD2)
axislegend(ax_t_MLD2, position = :rt)

ax_t_eR3 = Axis(fig[2, 3]; ylabel = "Tₑ (F₁₀₀₀/F₂₀₀)", xlabel = "t (day)", title = "monthly: Tₑ vs. NCP",
            yticklabelcolor=:darkgoldenrod2, ylabelcolor=:darkgoldenrod2)
lines!(ax_t_eR3, 1:1:365, vec(True_Teff2), color=:darkgoldenrod2, linewidth = 2, label = "Tₑ")
axislegend(ax_t_eR3, position = :lt)
ylims!(ax_t_eR3, 0, 1.3)
ax_t_MLD3 = Axis(fig[2, 3]; ylabel = "NCP (mmol m⁻³ d⁻¹)", xlabel = "t (day)", yaxisposition=:right, 
                yticklabelcolor=:royalblue1, ylabelcolor=:royalblue1)
lines!(ax_t_MLD3, 1:1:365, vec(tot_NCP2), linewidth = 2, color=:royalblue1, label = "NCP")
ylims!(ax_t_MLD3, 2,12)
linkxaxes!(ax_t_eR3, ax_t_MLD3)
axislegend(ax_t_MLD3, position = :rt)

display(fig)

#
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
# fig = Figure(size=(800, 500))

# ax_t_er = Axis(fig[1, 1]; ylabel = "ratio", xlabel = "t (day)", title = "e-ratio (F₂₀₀/NCP)")
# lines!(ax_t_er, 1:1:365, vec(True_eratio1), linewidth = 2, label = "Year")
# lines!(ax_t_er, 1:1:365, vec(True_eratio2), linewidth = 2, label = "Year+Month")
# # lines!(ax_t_er, 1:1:365, vec(True_eratio3), linewidth = 2, label = "Year+Week")
# ylims!(ax_t_er, 0,0.6)




