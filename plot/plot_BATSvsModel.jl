using GLMakie
using Printf
using Statistics
using XLSX

using Oceananigans
using Oceananigans.Units

xlsx_data = XLSX.readxlsx("bats_flux.xlsx")
data200 = xlsx_data["Sheet2"]
data300 = xlsx_data["Sheet3"]

# Extract the data from columns B and D
dayofyear_z200 = data200[2:end, 2]  # Column B (day)
POPflux_z200 = data200[2:end, 4]  # Column D (flux data)
sorted_ind200 = sortperm(vec(dayofyear_z200))
sorted_day_z200= dayofyear_z200[sorted_ind200]
sorted_flux_z200 = POPflux_z200[sorted_ind200]

dayofyear_z300 = data300[2:end, 2]  # Column B (day)
POPflux_z300 = data300[2:end, 4]
sorted_ind300 = sortperm(vec(dayofyear_z300))
sorted_day_z300= dayofyear_z300[sorted_ind300]
sorted_flux_z300 = POPflux_z300[sorted_ind300]

# bin the data into 10 days
bin_ind200 = ceil.(Int, sorted_day_z200 ./ 10)
unique_bins2 = unique(bin_ind200)
function round_to_nearest_5(x)
    return ceil(Int, x / 5) * 5
end
avg_t2 = [round_to_nearest_5.(minimum(sorted_day_z200[bin_ind200 .== bin])) for bin in unique_bins2]
avg_flux2 = [median(sorted_flux_z200[bin_ind200 .== bin]) for bin in unique_bins2]

bin_ind300 = ceil.(Int, sorted_day_z300 ./ 10)
unique_bins3 = unique(bin_ind300)
avg_t3 = [round_to_nearest_5.(minimum(sorted_day_z300[bin_ind300 .== bin])) for bin in unique_bins3]
avg_flux3 = [median(sorted_flux_z300[bin_ind300 .== bin]) for bin in unique_bins3]

# Load model results
filepath1 = "./P1_21y.jld2"
filepath2 = "./P2_21y.jld2"
filepath3 = "./P3_21y.jld2"
POP_timeseries1 = FieldTimeSeries(filepath1, "POP")
POP_timeseries2 = FieldTimeSeries(filepath2, "POP")
POP_timeseries3 = FieldTimeSeries(filepath3, "POP")
times = POP_timeseries1.times
xw, yw, zw = nodes(POP_timeseries1)

Oligo_2flux1 = zeros(365, 1) 
Oligo_2flux2 = zeros(365, 1) 
Oligo_2flux3 = zeros(365, 1)  
Oligo_3flux1 = zeros(365, 1) 
Oligo_3flux2 = zeros(365, 1) 
Oligo_3flux3 = zeros(365, 1)  

for i in 1:365
    # Compile 2D POP fluxes on each day
    Oligo_2flux1[i] = 1e4.*interior(POP_timeseries1[i], 1, 497, 190)
    Oligo_2flux2[i] = 1e4.*interior(POP_timeseries2[i], 1, 497, 190)
    Oligo_2flux3[i] = 1e4.*interior(POP_timeseries3[i], 1, 497, 190)
    Oligo_3flux1[i] = 1e4.*interior(POP_timeseries1[i], 1, 497, 185)
    Oligo_3flux2[i] = 1e4.*interior(POP_timeseries2[i], 1, 497, 185)
    Oligo_3flux3[i] = 1e4.*interior(POP_timeseries3[i], 1, 497, 185)
end

####################### Plot #######################
figure = Figure(size=(600, 1000))
axis = Axis(figure[1, 1], title="POP flux at 200 m", xlabel="t (day)", ylabel="POP flux (mmol m⁻² d⁻¹)")

lines!(axis, 1:1:365, vec(Oligo_2flux1), linewidth = 2, label = "Model - low")
lines!(axis, 1:1:365, vec(Oligo_2flux2), linewidth = 2, label = "Model - mid")
lines!(axis, 1:1:365, vec(Oligo_2flux3), linewidth = 2, label = "Model - high")

lines!(axis, avg_t2, avg_flux2, label = "BATS data")
ylims!(axis, 0, 0.01)
axislegend(axis, position = :rt)

axis2 = Axis(figure[2, 1], title="POP flux at 300 m", xlabel="t (day)", ylabel="POP flux (mmol m⁻² d⁻¹)")

lines!(axis2, 1:1:365, vec(Oligo_3flux1), linewidth = 2, label = "Model - low")
lines!(axis2, 1:1:365, vec(Oligo_3flux2), linewidth = 2, label = "Model - mid")
lines!(axis2, 1:1:365, vec(Oligo_3flux3), linewidth = 2, label = "Model - high")

lines!(axis2, avg_t3, avg_flux3, label = "BATS data")
ylims!(axis2, 0, 0.01)
axislegend(axis2, position = :rt)

display(figure)