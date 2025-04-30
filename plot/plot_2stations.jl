# Compare two stations in the 2D domain: productive vs oligotrophic

using GLMakie
using Printf
using Statistics

using Oceananigans
using Oceananigans.Units

filepath1 = "./P5_51y.jld2"
filepath2 = "./P6_51y.jld2"
# filepath3 = "./P3_21y.jld2"

NCP_timeseries1 = FieldTimeSeries(filepath1, "NCP")
NCP_timeseries2 = FieldTimeSeries(filepath2, "NCP")
# NCP_timeseries3 = FieldTimeSeries(filepath3, "NCP")

Premin_timeseries1 = FieldTimeSeries(filepath1, "Premin")
times = Premin_timeseries1.times
xw, yw, zw = nodes(Premin_timeseries1)
Premin_timeseries2 = FieldTimeSeries(filepath2, "Premin")
# Premin_timeseries3 = FieldTimeSeries(filepath3, "Premin")

POP_timeseries1 = FieldTimeSeries(filepath1, "POP")
POP_timeseries2 = FieldTimeSeries(filepath2, "POP")
# POP_timeseries3 = FieldTimeSeries(filepath3, "POP")

#################################################################
####### Timeseries POP flux: productive vs oligotrophic #######
#################################################################
dz = 20 

Prod_POP_flux1 = zeros(365, 200) 
Prod_POP_flux2 = zeros(365, 200) 
# Prod_POP_flux3 = zeros(365, 200)  
Oligo_POP_flux1 = zeros(365, 200) 
Oligo_POP_flux2 = zeros(365, 200) 
# Oligo_POP_flux3 = zeros(365, 200)  

Prod_eratio1 = zeros(365, 1) 
Prod_eratio2 = zeros(365, 1) 
# Prod_eratio3 = zeros(365, 1) 
Oligo_eratio1 = zeros(365, 1) 
Oligo_eratio2 = zeros(365, 1) 
# Oligo_eratio3 = zeros(365, 1) 

Prod_Teff1 = zeros(365, 1) 
Prod_Teff2 = zeros(365, 1) 
# Prod_Teff3 = zeros(365, 1) 
Oligo_Teff1 = zeros(365, 1) 
Oligo_Teff2 = zeros(365, 1) 
# Oligo_Teff3 = zeros(365, 1) 

for i in 1:365
    # Compile 2D POP fluxes on each day
    Prod_POP_flux1[i, :] = 1e4*interior(POP_timeseries1[i], 1, 115, :)
    Prod_POP_flux2[i, :] = 1e4*interior(POP_timeseries2[i], 1, 115, :)
    # Prod_POP_flux3[i, :] = 1e4*interior(POP_timeseries3[i], 1, 115, :)

    Oligo_POP_flux1[i, :] = 1e4*interior(POP_timeseries1[i], 1, 495, :)
    Oligo_POP_flux2[i, :] = 1e4*interior(POP_timeseries2[i], 1, 495, :)
    # Oligo_POP_flux3[i, :] = 1e4*interior(POP_timeseries3[i], 1, 495, :)

    # Calculate e-ratio (115 m/NCP)
    Prod_eratio1[i, 1] = Prod_POP_flux1[i, 195]/sum(dz*1days*1e3*interior(NCP_timeseries1[i], 1, 115, :))
    Prod_eratio2[i, 1] = Prod_POP_flux2[i, 195]/sum(dz*1days*1e3*interior(NCP_timeseries2[i], 1, 115, :))
    # Prod_eratio3[i, 1] = Prod_POP_flux3[i, 195]/sum(dz*1days*1e3*interior(NCP_timeseries3[i], 1, 115, :))
   
    Oligo_eratio1[i, 1] = Oligo_POP_flux1[i, 195]/sum(dz*1days*1e3*interior(NCP_timeseries1[i], 1, 495, :))
    Oligo_eratio2[i, 1] = Oligo_POP_flux2[i, 195]/sum(dz*1days*1e3*interior(NCP_timeseries2[i], 1, 495, :))
    # Oligo_eratio3[i, 1] = Oligo_POP_flux3[i, 195]/sum(dz*1days*1e3*interior(NCP_timeseries3[i], 1, 495, :))
   
    # Calculate transfer efficiency (1000 m/200 m)
    Prod_Teff1[i, 1] = sum(Prod_POP_flux1[i, 150])/sum(Prod_POP_flux1[i, 190])
    Prod_Teff2[i, 1] = sum(Prod_POP_flux2[i, 150])/sum(Prod_POP_flux2[i, 190])
    # Prod_Teff3[i, 1] = sum(Prod_POP_flux3[i, 150])/sum(Prod_POP_flux3[i, 190])

    Oligo_Teff1[i, 1] = sum(Oligo_POP_flux1[i, 150])/sum(Oligo_POP_flux1[i, 190])
    Oligo_Teff2[i, 1] = sum(Oligo_POP_flux2[i, 150])/sum(Oligo_POP_flux2[i, 190])
    # Oligo_Teff3[i, 1] = sum(Oligo_POP_flux3[i, 150])/sum(Oligo_POP_flux3[i, 190])
end

# Martin transfer efficiency (1000 m/200 m)
z₀ = log(0.01)*25 
Martin_value = ((zw[150]+z₀)/(zw[190]+z₀))^-0.84
Martin_Teff = repeat([Martin_value], 365, 1)

# t = 1:1:365
# MLD1 = 100 .+ 50 .* sinpi.(2 .* t ./ 365)
# MLD2 = 100 .+ 50 .* sinpi.(2 .* t ./ (365/4))
# MLD3 = 100 .+ 50 .* sinpi.(2 .* t ./ (365/12))

#################################################################
######################### Plot Teff together #########################
#################################################################

fig = Figure(size=(600, 800))

# ax_t_MLD = Axis(fig[1, 1]; ylabel = "MLD (m)", xlabel = "t (day)", 
#                 yticks=0:50:200)
# lines!(ax_t_MLD, 1:1:365, MLD1, linewidth = 2, label = "Low frequency")
# lines!(ax_t_MLD, 1:1:365, MLD2, linewidth = 2, label = "Mid frequency")
# lines!(ax_t_MLD, 1:1:365, MLD3, linewidth = 2, label = "High frequency")
# ylims!(ax_t_MLD, 200,0)
# axislegend(ax_t_MLD, position = :rt)
# display(fig)
#
ax_t_prod = Axis(fig[1, 2]; ylabel = "Ratio of F₁₀₀₀/F₂₀₀", xlabel = "t (day)", title = "Transfer effeciency (productive)")
lines!(ax_t_prod, 1:1:365, vec(Prod_Teff1), linewidth = 2, label = "Year")
lines!(ax_t_prod, 1:1:365, vec(Prod_Teff2), linewidth = 2, label = "Month")
# lines!(ax_t_prod, 1:1:365, vec(Prod_Teff3), linewidth = 2, label = "Year+Week")
lines!(ax_t_prod, 1:1:365,vec(Martin_Teff),linewidth = 2, label = "Martin")
axislegend(ax_t_prod, position = :lt)
ylims!(ax_t_prod, 0, 1)

ax_t_oligo = Axis(fig[2, 2]; ylabel = "Ratio of F₁₀₀₀/F₂₀₀", xlabel = "t (day)", title = "Transfer effeciency (oligotrophic)")
lines!(ax_t_oligo, 1:1:365, vec(Oligo_Teff1), linewidth = 2, label = "Year")
lines!(ax_t_oligo, 1:1:365, vec(Oligo_Teff2), linewidth = 2, label = "Month")
# lines!(ax_t_oligo, 1:1:365, vec(Oligo_Teff3), linewidth = 2, label = "Year+Week")
lines!(ax_t_oligo, 1:1:365,vec(Martin_Teff),linewidth = 2, label = "Martin")
# axislegend(ax_t_oligo, position = :rb)
ylims!(ax_t_oligo, 0, 1)
# display(fig)
# export ratio
ax_r_prod = Axis(fig[1, 1]; ylabel = "F₁₁₅/NCP", xlabel = "t (day)", title = "e-ratio (productive)")
lines!(ax_r_prod, 1:1:365, vec(Prod_eratio1), linewidth = 2, label = "Year")
lines!(ax_r_prod, 1:1:365, vec(Prod_eratio2), linewidth = 2, label = "Month")
# lines!(ax_r_prod, 1:1:365, vec(Prod_eratio3), linewidth = 2, label = "Model - high")
# axislegend(ax_r_prod, position = :rb)
ylims!(ax_r_prod, 0, 0.5)

ax_r_oligo = Axis(fig[2, 1]; ylabel = "F₁₁₅/NCP", xlabel = "t (day)", title = "e-ratio (oligotrophic)")
lines!(ax_r_oligo, 1:1:365, vec(Oligo_eratio1), linewidth = 2, label = "Year")
lines!(ax_r_oligo, 1:1:365, vec(Oligo_eratio2), linewidth = 2, label = "Month")
# lines!(ax_r_oligo, 1:1:365, vec(Oligo_eratio3), linewidth = 2, label = "Model - high")
# axislegend(ax_r_oligo, position = :rb)
ylims!(ax_r_oligo, 0, 0.5)

display(fig)
