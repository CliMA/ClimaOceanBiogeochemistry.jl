using GLMakie
using Printf
using Statistics

using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients
using Oceananigans
using Oceananigans.Units
using Oceananigans.Fields: ZeroField, CenterField
using Oceananigans.BoundaryConditions: fill_halo_regions!

filepath1 = "./1perY_11y.jld2"
filepath2 = "./4perY_11y.jld2"
filepath3 = "./12perY_11y.jld2"

#################################################################
##### Load data: rates of NCP, Dremin, Premin #####
#################################################################

NCP_timeseries1 = FieldTimeSeries(filepath1, "NCP")
NCP_timeseries2 = FieldTimeSeries(filepath2, "NCP")
NCP_timeseries3 = FieldTimeSeries(filepath3, "NCP")

POP_timeseries1 = FieldTimeSeries(filepath1, "POP")
POP_timeseries2 = FieldTimeSeries(filepath2, "POP")
POP_timeseries3 = FieldTimeSeries(filepath3, "POP")

times = POP_timeseries1.times
xw, yw, zw = nodes(POP_timeseries1)

#################################################################
######### prod vs. remin #########
#################################################################
#################################################################
####### Timeseries POP flux: productive vs oligotrophic #######
#################################################################
dz = 20 

Prod_NCP1 = zeros(365, 1) 
Prod_NCP2 = zeros(365, 1) 
Prod_NCP3 = zeros(365, 1)  

Oligo_NCP1 = zeros(365, 1) 
Oligo_NCP2 = zeros(365, 1) 
Oligo_NCP3 = zeros(365, 1)  

Prod_POP_flux1 = zeros(365, 200) 
Prod_POP_flux2 = zeros(365, 200) 
Prod_POP_flux3 = zeros(365, 200)  

Oligo_POP_flux1 = zeros(365, 200) 
Oligo_POP_flux2 = zeros(365, 200) 
Oligo_POP_flux3 = zeros(365, 200)  

Prod_100flux1 = zeros(365, 1) 
Prod_100flux2 = zeros(365, 1) 
Prod_100flux3 = zeros(365, 1)  
Oligo_100flux1 = zeros(365, 1) 
Oligo_100flux2 = zeros(365, 1) 
Oligo_100flux3 = zeros(365, 1)   

Prod_1000flux1 = zeros(365, 1) 
Prod_1000flux2 = zeros(365, 1) 
Prod_1000flux3 = zeros(365, 1)  
Oligo_1000flux1 = zeros(365, 1) 
Oligo_1000flux2 = zeros(365, 1) 
Oligo_1000flux3 = zeros(365, 1)  

Prod_1000martin1 = zeros(365, 1) 
Prod_1000martin2 = zeros(365, 1) 
Prod_1000martin3 = zeros(365, 1)  
Oligo_1000martin1 = zeros(365, 1) 
Oligo_1000martin2 = zeros(365, 1) 
Oligo_1000martin3 = zeros(365, 1)  

z₀ = log(0.01)*25 
Martin_factor = ((zw[150].+z₀)./(zw[190]+z₀)).^-0.84

for i in 1:365
    # Calculate NCP (mmol m⁻³ d⁻¹)
    Prod_NCP1[i, 1] = sum(1days*1e3*interior(NCP_timeseries1[i], 1, 115, :))
    Prod_NCP2[i, 1] = sum(1days*1e3*interior(NCP_timeseries2[i], 1, 115, :))
    Prod_NCP3[i, 1] = sum(1days*1e3*interior(NCP_timeseries3[i], 1, 115, :))
    
    Oligo_NCP1[i, 1] = sum(1days*1e3*interior(NCP_timeseries1[i], 1, 490, :))
    Oligo_NCP2[i, 1] = sum(1days*1e3*interior(NCP_timeseries2[i], 1, 490, :))
    Oligo_NCP3[i, 1] = sum(1days*1e3*interior(NCP_timeseries3[i], 1, 490, :))

    # Compile 2D POP fluxes on each day
    Prod_POP_flux1[i, :] = 1e4*interior(POP_timeseries1[i], 1, 115, :)
    Prod_POP_flux2[i, :] = 1e4*interior(POP_timeseries2[i], 1, 115, :)
    Prod_POP_flux3[i, :] = 1e4*interior(POP_timeseries3[i], 1, 115, :)

    Oligo_POP_flux1[i, :] = 1e4*interior(POP_timeseries1[i], 1, 490, :)
    Oligo_POP_flux2[i, :] = 1e4*interior(POP_timeseries2[i], 1, 490, :)
    Oligo_POP_flux3[i, :] = 1e4*interior(POP_timeseries3[i], 1, 490, :)

    # POP flux at 100 m
    Prod_100flux1[i] = Prod_POP_flux1[i, 195] 
    Prod_100flux2[i] = Prod_POP_flux2[i, 195] 
    Prod_100flux3[i] = Prod_POP_flux3[i, 195] 

    Oligo_100flux1[i] = Oligo_POP_flux1[i, 195] 
    Oligo_100flux2[i] = Oligo_POP_flux2[i, 195] 
    Oligo_100flux3[i] = Oligo_POP_flux3[i, 195] 

    # POP flux at 1000 m
    Prod_1000flux1[i] = Prod_POP_flux1[i, 150] 
    Prod_1000flux2[i] = Prod_POP_flux2[i, 150] 
    Prod_1000flux3[i] = Prod_POP_flux3[i, 150] 

    Oligo_1000flux1[i] = Oligo_POP_flux1[i, 150] 
    Oligo_1000flux2[i] = Oligo_POP_flux2[i, 150] 
    Oligo_1000flux3[i] = Oligo_POP_flux3[i, 150] 

    # MArtin-POP flux at 1000 m
    Prod_1000martin1[i] = Prod_POP_flux1[i, 190] * Martin_factor 
    Prod_1000martin2[i] = Prod_POP_flux2[i, 190] * Martin_factor 
    Prod_1000martin3[i] = Prod_POP_flux3[i, 190] * Martin_factor 
    
    Oligo_1000martin1[i] = Oligo_POP_flux1[i, 190] * Martin_factor  
    Oligo_1000martin2[i] = Oligo_POP_flux2[i, 190] * Martin_factor  
    Oligo_1000martin3[i] = Oligo_POP_flux3[i, 190] * Martin_factor  

end

#################################################################
######################### Plot  #########################
#################################################################

fig = Figure(size=(1500, 1200))

ax_ncp_prod = Axis(fig[1, 1]; ylabel = "NCP (mmol m⁻³ d⁻¹)", xlabel = "t (day)", title = "NCP (productive)")
lines!(ax_ncp_prod, 1:1:365, vec(Prod_NCP1), linewidth = 2, label = "low")
lines!(ax_ncp_prod, 1:1:365, vec(Prod_NCP2), linewidth = 2, label = "mid")
lines!(ax_ncp_prod, 1:1:365, vec(Prod_NCP3), linewidth = 2, label = "high")
axislegend(ax_ncp_prod, position = :rb)
ylims!(ax_ncp_prod, 0.0185, 0.0235)

ax_ncp_oligo = Axis(fig[1, 2]; ylabel = "NCP (mmol m⁻³ d⁻¹)", xlabel = "t (day)", title = "NCP (oligotrophic)")
lines!(ax_ncp_oligo, 1:1:365, vec(Oligo_NCP1), linewidth = 2, label = "low")
lines!(ax_ncp_oligo, 1:1:365, vec(Oligo_NCP2), linewidth = 2, label = "mid")
lines!(ax_ncp_oligo, 1:1:365, vec(Oligo_NCP3), linewidth = 2, label = "high")
axislegend(ax_ncp_oligo, position = :rb)
ylims!(ax_ncp_oligo, 1.76e-3, 1.83e-3)

ax_100_prod = Axis(fig[2, 1]; ylabel = "POP flux (mmol m⁻² d⁻¹)", xlabel = "t (day)", title = "POP flux at 100 m (productive)")
lines!(ax_100_prod, 1:1:365, vec(Prod_100flux1), linewidth = 2, label = "low")
lines!(ax_100_prod, 1:1:365, vec(Prod_100flux2), linewidth = 2, label = "mid")
lines!(ax_100_prod, 1:1:365, vec(Prod_100flux3), linewidth = 2, label = "high")
axislegend(ax_100_prod, position = :rb)
ylims!(ax_100_prod, 0.075, 0.095)

ax_100_oligo = Axis(fig[2, 2]; ylabel = "POP flux (mmol m⁻² d⁻¹)", xlabel = "t (day)", title = "POP flux at 100 m (oligotrophic)")
lines!(ax_100_oligo, 1:1:365, vec(Oligo_100flux1), linewidth = 2, label = "low")
lines!(ax_100_oligo, 1:1:365, vec(Oligo_100flux2), linewidth = 2, label = "mid")
lines!(ax_100_oligo, 1:1:365, vec(Oligo_100flux3), linewidth = 2, label = "high")
axislegend(ax_100_oligo, position = :rb)
ylims!(ax_100_oligo, 0.0069, 0.0076)

ax_1000_prod = Axis(fig[3, 1]; ylabel = "POP flux (mmol m⁻² d⁻¹)", xlabel = "t (day)", title = "POP flux at 1000 m (productive)")
lines!(ax_1000_prod, 1:1:365, vec(Prod_1000flux1), linewidth = 2, label = "low")
lines!(ax_1000_prod, 1:1:365, vec(Prod_1000flux2), linewidth = 2, label = "mid")
lines!(ax_1000_prod, 1:1:365, vec(Prod_1000flux3), linewidth = 2, label = "high")
axislegend(ax_1000_prod, position = :rb)
ylims!(ax_1000_prod, 0.022, 0.028)

ax_1000_oligo = Axis(fig[3, 2]; ylabel = "POP flux (mmol m⁻² d⁻¹)", xlabel = "t (day)", title = "POP flux at 1000 m (oligotrophic)")
lines!(ax_1000_oligo, 1:1:365, vec(Oligo_1000flux1), linewidth = 2, label = "low")
lines!(ax_1000_oligo, 1:1:365, vec(Oligo_1000flux2), linewidth = 2, label = "mid")
lines!(ax_1000_oligo, 1:1:365, vec(Oligo_1000flux3), linewidth = 2, label = "high")
axislegend(ax_1000_oligo, position = :rb)
ylims!(ax_1000_oligo, 0.00195, 0.00219)

ax_1000m_prod = Axis(fig[4, 1]; ylabel = "POP flux (mmol m⁻² d⁻¹)", xlabel = "t (day)", title = "Martin-calculated flux at 1000 m (productive)")
lines!(ax_1000m_prod, 1:1:365, vec(Prod_1000martin1), linewidth = 2, label = "low")
lines!(ax_1000m_prod, 1:1:365, vec(Prod_1000martin2), linewidth = 2, label = "mid")
lines!(ax_1000m_prod, 1:1:365, vec(Prod_1000martin3), linewidth = 2, label = "high")
axislegend(ax_1000m_prod, position = :rb)
ylims!(ax_1000m_prod, 0.022, 0.028)

ax_1000m_oligo = Axis(fig[4, 2]; ylabel = "POP flux (mmol m⁻² d⁻¹)", xlabel = "t (day)", title = "Martin-calculated flux at 1000 m (oligotrophic)")
lines!(ax_1000m_oligo, 1:1:365, vec(Oligo_1000martin1), linewidth = 2, label = "low")
lines!(ax_1000m_oligo, 1:1:365, vec(Oligo_1000martin2), linewidth = 2, label = "mid")
lines!(ax_1000m_oligo, 1:1:365, vec(Oligo_1000martin3), linewidth = 2, label = "high")
axislegend(ax_1000m_oligo, position = :rb)
ylims!(ax_1000m_oligo, 0.00195, 0.00219)

display(fig)


####################### flux vs Martin ####################### 
fig2 = Figure(size=(1200, 1200))

ax_1_prod = Axis(fig2[1, 1]; ylabel = "POP flux (mmol m⁻² d⁻¹)", xlabel = "t (day)", title = "Low frequency (productive)")
lines!(ax_1_prod, 1:1:365, vec(Prod_1000flux1), linewidth = 2, label = "model")
lines!(ax_1_prod, 1:1:365, vec(Prod_1000martin1), linewidth = 2, label = "Martin")
axislegend(ax_1_prod, position = :rb)
ylims!(ax_1_prod, 0.022, 0.028)

ax_1_oligo = Axis(fig2[1, 2]; ylabel = "POP flux (mmol m⁻² d⁻¹)", xlabel = "t (day)", title = "Low frequency (oligotrophic)")
lines!(ax_1_oligo, 1:1:365, vec(Oligo_1000flux1), linewidth = 2, label = "model")
lines!(ax_1_oligo, 1:1:365, vec(Oligo_1000martin1), linewidth = 2, label = "Martin")
axislegend(ax_1_oligo, position = :rb)
ylims!(ax_1_oligo, 0.0019, 0.0022)

ax_2_prod = Axis(fig2[2, 1]; ylabel = "POP flux (mmol m⁻² d⁻¹)", xlabel = "t (day)", title = "Mid frequency (productive)")
lines!(ax_2_prod, 1:1:365, vec(Prod_1000flux2), linewidth = 2, label = "model")
lines!(ax_2_prod, 1:1:365, vec(Prod_1000martin2), linewidth = 2, label = "Martin")
axislegend(ax_2_prod, position = :rb)
ylims!(ax_2_prod, 0.022, 0.028)

ax_2_oligo = Axis(fig2[2, 2]; ylabel = "POP flux (mmol m⁻² d⁻¹)", xlabel = "t (day)", title = "Mid frequency (oligotrophic)")
lines!(ax_2_oligo, 1:1:365, vec(Oligo_1000flux2), linewidth = 2, label = "model")
lines!(ax_2_oligo, 1:1:365, vec(Oligo_1000martin2), linewidth = 2, label = "Martin")
axislegend(ax_2_oligo, position = :rb)
ylims!(ax_2_oligo, 0.0019, 0.0022)

ax_3_prod = Axis(fig2[3, 1]; ylabel = "POP flux (mmol m⁻² d⁻¹)", xlabel = "t (day)", title = "High frequency (productive)")
lines!(ax_3_prod, 1:1:365, vec(Prod_1000flux3), linewidth = 2, label = "model")
lines!(ax_3_prod, 1:1:365, vec(Prod_1000martin3), linewidth = 2, label = "Martin")
axislegend(ax_3_prod, position = :rb)
ylims!(ax_3_prod, 0.022, 0.028)

ax_3_oligo = Axis(fig2[3, 2]; ylabel = "POP flux (mmol m⁻² d⁻¹)", xlabel = "t (day)", title = "High frequency (oligotrophic)")
lines!(ax_3_oligo, 1:1:365, vec(Oligo_1000flux3), linewidth = 2, label = "model")
lines!(ax_3_oligo, 1:1:365, vec(Oligo_1000martin3), linewidth = 2, label = "Martin")
axislegend(ax_3_oligo, position = :rb)
ylims!(ax_3_oligo, 0.0019, 0.0022)

display(fig2)