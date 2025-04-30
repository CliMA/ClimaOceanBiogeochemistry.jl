using GLMakie
using Printf
using Statistics
using XLSX

using Oceananigans
using Oceananigans.Units

################################################
########### BATS Primary Production ############
################################################
xlsx_data = XLSX.readxlsx("bats_flux.xlsx")
OC_150 = xlsx_data["C_avg_150m"]
OC_200 = xlsx_data["C_avg_200m"]
OC_300 = xlsx_data["C_avg_300m"]

# primary production
pp_data = XLSX.readxlsx("bats_primary_production.xlsx")
data_pp = pp_data["Sheet1"]

dayofyear_pp = data_pp[2:end, 2]  # Column B (day)
bats_pp = data_pp[2:end, 3]  # Column C (flux data)
sorted_indpp = sortperm(vec(dayofyear_pp))
sorted_day_pp= dayofyear_pp[sorted_indpp]
sorted_pp = bats_pp[sorted_indpp]

unique_day_pp = unique(sorted_day_pp)
pp_mean = Float64[]
pp_stderr = Float64[]
for d in unique_day_pp
    idx = findall(x -> x == d, sorted_day_pp)
    values = sorted_pp[idx]
    push!(pp_mean, mean(values))
    push!(pp_stderr, std(values) / sqrt(length(values)))  # standard error
end

# Extract the data from columns B and D
dayofyear_z150 = OC_150[2:end, 3]  
OCflux_z150 = OC_150[2:end, 4]  
sorted_ind150 = sortperm(vec(dayofyear_z150))
sorted_day_z150= dayofyear_z150[sorted_ind150]
sorted_flux_z150 = OCflux_z150[sorted_ind150]

# There is one outlier
idx_outlier = findfirst(>(800), sorted_flux_z150)
sorted_day_z150 = deleteat!(copy(sorted_day_z150), idx_outlier)
sorted_flux_z150 = deleteat!(copy(sorted_flux_z150), idx_outlier)

# Calculate daily mean
unique_day_f150 = unique(sorted_day_z150)
f150_mean = Float64[]
f150_stderr = Float64[]
for d in unique_day_f150
    idx = findall(x -> x == d, sorted_day_z150)
    values = sorted_flux_z150[idx]
    push!(f150_mean, mean(values))
    push!(f150_stderr, std(values) / sqrt(length(values)))  # standard error
end

dayofyear_z200 = OC_200[2:end, 3]  
OCflux_z200 = OC_200[2:end, 4]  
sorted_ind200 = sortperm(vec(dayofyear_z200))
sorted_day_z200= dayofyear_z200[sorted_ind200]
sorted_flux_z200 = OCflux_z200[sorted_ind200]
# There is one outlier
idx_outlier = findall(>(300), sorted_flux_z200)
sorted_day_z200 = deleteat!(copy(sorted_day_z200), idx_outlier)
sorted_flux_z200 = deleteat!(copy(sorted_flux_z200), idx_outlier)
# Calculate daily mean
unique_day_f200 = unique(sorted_day_z200)
f200_mean = Float64[]
f200_stderr = Float64[]
for d in unique_day_f200
    idx = findall(x -> x == d, sorted_day_z200)
    values = sorted_flux_z200[idx]
    push!(f200_mean, mean(values))
    push!(f200_stderr, std(values) / sqrt(length(values)))  # standard error
end

dayofyear_z300 = OC_300[2:end, 3]  
OCflux_z300 = OC_300[2:end, 4]
sorted_ind300 = sortperm(vec(dayofyear_z300))
sorted_day_z300= dayofyear_z300[sorted_ind300]
sorted_flux_z300 = OCflux_z300[sorted_ind300]
# There is one outlier
idx_outlier = findfirst(>(400), sorted_flux_z300)
sorted_day_z300 = deleteat!(copy(sorted_day_z300), idx_outlier)
sorted_flux_z300 = deleteat!(copy(sorted_flux_z300), idx_outlier)
# Calculate daily mean
unique_day_f300 = unique(sorted_day_z300)
f300_mean = Float64[]
f300_stderr = Float64[]
for d in unique_day_f300
    idx = findall(x -> x == d, sorted_day_z300)
    values = sorted_flux_z300[idx]
    push!(f300_mean, mean(values))
    push!(f300_stderr, std(values) / sqrt(length(values)))  # standard error
end

#= bin the data into 10 days
function round_to_nearest_5(x)
    return ceil(Int, x / 5) * 5
end

bin_indpp = ceil.(Int, sorted_day_pp ./ 10)
unique_bins = unique(bin_indpp)
avg_t = [round_to_nearest_5.(minimum(sorted_day_pp[bin_indpp .== bin])) for bin in unique_bins]
avg_pp = [median(sorted_pp[bin_indpp .== bin]) for bin in unique_bins]
# Unit conversion: mg C m-3 d-1 to mmol P m-3 d-1
# ave_pp_P = avg_pp ./ 12 ./117

bin_ind200 = ceil.(Int, sorted_day_z200 ./ 10)
unique_bins2 = unique(bin_ind200)
avg_t2 = [round_to_nearest_5.(minimum(sorted_day_z200[bin_ind200 .== bin])) for bin in unique_bins2]
avg_flux2 = [median(sorted_flux_z200[bin_ind200 .== bin]) for bin in unique_bins2]

bin_ind300 = ceil.(Int, sorted_day_z300 ./ 10)
unique_bins3 = unique(bin_ind300)
avg_t3 = [round_to_nearest_5.(minimum(sorted_day_z300[bin_ind300 .== bin])) for bin in unique_bins3]
avg_flux3 = [median(sorted_flux_z300[bin_ind300 .== bin]) for bin in unique_bins3]
=#
#
figure = Figure(size=(1000, 600))

### Primary production ###
# ax_pp1 = Axis(figure[1, 1], title="Primary production", xlabel="t (day)", ylabel="Productivity (mg C m⁻³ d⁻¹)")
# scatter!(ax_pp1,sorted_day_pp, sorted_pp; color=:grey18, alpha = 0.4, markersize = 8)

ax_pp2 = Axis(figure[1, 1], title="(a) Primary production", xlabel="t (day)", ylabel="Productivity (mg C m⁻³ d⁻¹)")
band!(ax_pp2, 0.0:1:365.0, fill(mean(pp_mean).-(std(pp_mean)), 366), fill(mean(pp_mean).+(std(pp_mean)), 366); color = (:gray, 0.2))
lines!(ax_pp2, unique_day_pp, pp_mean, color=:grey25)
# errorbars!(ax_pp2, Float64.(unique_day_pp), pp_mean, ylower=pp_stderr, yupper=pp_stderr, whiskerwidth=5, color=:blue)

# ax_f150 = Axis(figure[2, 1], title="OC flux at 150 m", xlabel="t (day)", ylabel="OC flux (mg m⁻² d⁻¹)")
# scatter!(ax_f150,sorted_day_z150, sorted_flux_z150; color=:grey18, alpha = 0.4, markersize = 8)

ax_f150_2 = Axis(figure[1, 2], title="(b) OC flux at 150 m", xlabel="t (day)", ylabel="OC flux (mg m⁻² d⁻¹)")
band!(ax_f150_2, 0.0:1:365.0, fill(mean(f150_mean).-(std(f150_mean)), 366), fill(mean(f150_mean).+(std(f150_mean)), 366); color = (:gray, 0.2))
lines!(ax_f150_2, unique_day_f150, f150_mean, color=:grey25)

#=
ax_f200 = Axis(figure[3, 1], title="OC flux at 200 m", xlabel="t (day)", ylabel="OC flux (mg m⁻² d⁻¹)")
scatter!(ax_f200,sorted_day_z200, sorted_flux_z200; color=:grey18, alpha = 0.4, markersize = 8)

ax_f200_2 = Axis(figure[3, 2], title="OC flux at 200 m (daily mean)", xlabel="t (day)", ylabel="OC flux (mg m⁻² d⁻¹)")
band!(ax_f200_2, 0.0:1:365.0, fill(mean(f200_mean)/2, 366), fill(mean(f200_mean)*2, 366); color = (:gray, 0.2))
lines!(ax_f200_2, unique_day_f200, f200_mean, color=:grey25)

ax_f300 = Axis(figure[4, 1], title="OC flux at 300 m", xlabel="t (day)", ylabel="OC flux (mg m⁻² d⁻¹)")
scatter!(ax_f300,sorted_day_z300, sorted_flux_z300; color=:grey18, alpha = 0.4, markersize = 8)

ax_f300_2 = Axis(figure[4, 2], title="OC flux at 300 m (daily mean)", xlabel="t (day)", ylabel="OC flux (mg m⁻² d⁻¹)")
band!(ax_f300_2, 0.0:1:365.0, fill(mean(f300_mean)/2, 366), fill(mean(f300_mean)*2, 366); color = (:gray, 0.2))
lines!(ax_f300_2, unique_day_f300, f300_mean, color=:grey25)

display(figure)
=#
########## Plot export ratio ############
common_times = intersect(unique_day_pp, unique_day_f150)
pp_common = [pp_mean[argmin(abs.(unique_day_pp .- t))] for t in common_times]
f150_common = [f150_mean[argmin(abs.(unique_day_f150 .- t))] for t in common_times]
e_ratio = f150_common ./ pp_common/150 # 150 for depth correction (unit)

relative_anomaly_eratio = (e_ratio .- mean(e_ratio)) ./ mean(e_ratio)

# figure2 = Figure(size=(600, 300))
ax_r = Axis(figure[1, 3], title="(c) Export ratio = F₁₅₀/PP", xlabel="t (day)", ylabel="e-ratio")
lines!(ax_r, common_times, e_ratio, color=:black)
band!(ax_r, 0.0:1:365.0, fill(mean(e_ratio).-std(e_ratio), 366), fill(mean(e_ratio).+std(e_ratio), 366); color = (:gray, 0.2))

# ax_a = Axis(figure[2, 2], title="% anomaly of F₁₅₀/PP", xlabel="t (day)", ylabel="e-ratio anomaly (%)")
# lines!(ax_a, common_times, relative_anomaly_eratio .*100, color=:black)

display(figure)
#
#=
### OC flux ###
axis = Axis(figure[2, 1], title="OC flux at 200 m", xlabel="t (day)", ylabel="OC flux (mg m⁻² d⁻¹)")
lines!(axis, avg_t2, avg_flux2, color=:grey43, label = "BATS POP flux")
ylims!(axis, 0, 0.013)

axis2 = Axis(figure[2, 2], title="OC flux at 300 m", xlabel="t (day)", ylabel="OC flux (mg m⁻² d⁻¹)")
lines!(axis2, avg_t3, avg_flux3, color=:grey43, label = "BATS POP flux")
ylims!(axis2, 0, 0.01)
# axislegend(axis2, position = :rt)

display(figure)
=#