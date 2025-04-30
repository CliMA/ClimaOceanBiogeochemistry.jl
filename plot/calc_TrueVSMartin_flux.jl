using GLMakie
using Printf
using Statistics

using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients
using Oceananigans
using Oceananigans.Units
using Oceananigans.Fields: ZeroField, CenterField
using Oceananigans.BoundaryConditions: fill_halo_regions!

filepath1 = "./P6_51y.jld2"
NCP_timeseries1 = FieldTimeSeries(filepath1, "NCP")
Premin_timeseries1 = FieldTimeSeries(filepath1, "Premin")
times = Premin_timeseries1.times
xw, yw, zw = nodes(Premin_timeseries1)
POP_timeseries1 = FieldTimeSeries(filepath1, "POP")


dz = 20 # z grid gap
z₀ = log(0.01)*25 
Martin_factor = ((zw[190] .+ z₀) ./ (zw .+ z₀)).^0.84

######################################################################
########################## calculate fluxes ##########################
######################################################################
#
POP_flux1 = zeros(365, 500, 200) 

for i in 1:365
    # Compile 2D POP fluxes on each day
    POP_flux1[i, :, :] = 1e4*interior(POP_timeseries1[i], 1, :, :)
end
sequester_Ftrue = sum(POP_flux1[:, :, 150]; dims=2) # flux at 1000 m
export_Ftrue = sum(POP_flux1[:, :, 175]; dims=2)
deep_Ftrue = sum(POP_flux1[:, :, 50]; dims=2) 

POP_ref1 = POP_flux1[:, :, 190]
# POP flux in Martin curve (mmol m-2 d-1)
POP_Martin1 = zeros(365, 500, 200) 
for i in 1:365
    temp_Martin1 = repeat(POP_ref1[i,:], 1, 200)
    for k in 1:200
        temp_Martin1[:,k] .= temp_Martin1[:,k] .* Martin_factor[k]
    end
    POP_Martin1[i, :, :] = temp_Martin1
end
sequester_Fmartin = sum(POP_Martin1[:, :, 150]; dims=2) 
export_Fmartin = sum(POP_Martin1[:, :, 175]; dims=2) 
deep_Fmartin = sum(POP_Martin1[:, :, 50]; dims=2) 

#################################################################
######################### Annual flux estimates #########################
#################################################################

fig = Figure(size=(1200, 500))

ax_t1 = Axis(fig[1, 1]; ylabel = "flux (mmol m⁻² d⁻¹)", xlabel = "t (day)", title = "Integrated flux at 500 m in a year")
lines!(ax_t1, 1:1:365, cumsum(vec(export_Ftrue)), linewidth = 2, label = "Model")
lines!(ax_t1, 1:1:365,cumsum(vec(export_Fmartin)),linewidth = 2, label = "Martin")
axislegend(ax_t1, position = :rb)

ax_t2 = Axis(fig[1, 2]; ylabel = "flux (mmol m⁻² d⁻¹)", xlabel = "t (day)", title = "Integrated flux at 1000 m in a year")
lines!(ax_t2, 1:1:365, cumsum(vec(sequester_Ftrue)), linewidth = 2, label = "Model")
lines!(ax_t2, 1:1:365,cumsum(vec(sequester_Fmartin)),linewidth = 2, label = "Martin")
axislegend(ax_t2, position = :rb)
# ylims!(ax_t2, 0.3, 0.4)

ax_t3 = Axis(fig[1, 3]; ylabel = "flux (mmol m⁻² d⁻¹)", xlabel = "t (day)", title = "Integrated flux at 3000 m in a year")
lines!(ax_t3, 1:1:365, cumsum(vec(deep_Ftrue)), linewidth = 2, label = "Model")
lines!(ax_t3, 1:1:365,cumsum(vec(deep_Fmartin)),linewidth = 2, label = "Martin")
axislegend(ax_t3, position = :rb)

display(fig)