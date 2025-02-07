using GLMakie
using Printf
using Statistics

using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients
using Oceananigans
using Oceananigans.Units
using Oceananigans.Fields: ZeroField, CenterField
using Oceananigans.BoundaryConditions: fill_halo_regions!

Ny = 500 
Nz = 200
Ly = 15000kilometers   # m
Lz = 4000           # m

arch = CPU()
# We use a two-dimensional grid, with a `Flat` `y`-direction:
grid = RectilinearGrid(arch,
                       size = (Ny, Nz),
                       y = (0, Ly),
                       z = (-Lz, 0),
                       topology=(Flat, Bounded, Bounded))

# Load model outputs
filepath = "./AMOC115_auxiliary.jld2"

# using JLD2
# data = load("AMOC115_auxiliary.jld2")
# println(keys(data))  # List variable names

PO4_timeseries = FieldTimeSeries(filepath, "PO₄")
xt, yt, zt = nodes(PO4_timeseries)
POP_timeseries = FieldTimeSeries(filepath, "POP")
DOP_timeseries = FieldTimeSeries(filepath, "DOP")
Fe_timeseries = FieldTimeSeries(filepath, "Fe")
NO3_timeseries = FieldTimeSeries(filepath, "NO₃")
NCP_timeseries = FieldTimeSeries(filepath, "NCP")
Premin_timeseries = FieldTimeSeries(filepath, "Premin")
Dremin_timeseries = FieldTimeSeries(filepath, "Dremin")

# Define color maps
colors5 = cgrad(:RdYlBu, 5, categorical=true, rev=true) 
colors6 = cgrad(:RdYlBu, 6, categorical=true, rev=true) 
colors8 = cgrad(:RdYlBu, 8, categorical=true, rev=true)

# Load final time point
# tracer concentration unit: mol m⁻³ 
PO4_final = interior(PO4_timeseries[end], 1, :, :)
NO3_final = interior(NO3_timeseries[end], 1, :, :)
POP_final = interior(POP_timeseries[end], 1, :, :)
DOP_final = interior(DOP_timeseries[end], 1, :, :)
Fe_final = interior(Fe_timeseries[end], 1, :, :)

NCP_final = interior(NCP_timeseries[end], 1, :, :) .* 1e3 .* 1days
Premin_final = interior(Premin_timeseries[end], 1, :, :) .* 1e3 .* 1days
Dremin_final = interior(Dremin_timeseries[end], 1, :, :) .* 1e3 .* 1days

# Load model parameters from CAN
z_matrix = repeat(zt, 1, 500)
kᴵ=30
kᴾ=5e-7 * 1024.5  # mol P m⁻³ 
kᴺ=7e-6 * 1024.5 
kᶠ=1e-10 * 1024.5  

# light field
incident_PAR =  Field{Nothing, Center, Center}(grid)
surface_PAR(y,z) = 700 * sinpi(y/Ly) 
set!(incident_PAR, surface_PAR)   
fill_halo_regions!(incident_PAR, arch)

I = incident_PAR[1,1:500,1:200] .* exp.(z_matrix' ./ 25)
P = max.(0, PO4_final)
N = max.(0, NO3_final)
F = max.(0, Fe_final)

# Limitation terms
light_lim = I ./ (I .+ kᴵ)
p_lim = P ./ (P .+ kᴾ)
n_lim = N ./ (N .+ kᴺ)
f_lim = F ./ (F .+ kᶠ)
all_lim = light_lim .* min.(p_lim, n_lim, f_lim)  

#= tracer concentration and limitation terms
fig_can = Figure(size = (1100, 600))

# Zoom into surface
ax_l_lim2 = Axis(fig_can[1, 1]; xlabel = "y (km)", ylabel = "z (m)",title = "Surface I/(I+kᵢ)", aspect = 1)
hm_l_lim2 = heatmap!(ax_l_lim2, yt/1e3, zt[190:200], light_lim[:,190:200]; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_can[1, 2], hm_l_lim2; flipaxis = false)
# NO3 limitation
ax_n_lim = Axis(fig_can[2, 1]; xlabel = "y (km)", ylabel = "z (m)",title = "NO₃/(NO₃+kₙₒ₃)", aspect = 1)
hm_n_lim = heatmap!(ax_n_lim, yt/1e3, zt, n_lim; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_can[2, 2], hm_n_lim; flipaxis = false)
# PO4 limitation
ax_p_lim = Axis(fig_can[2, 3]; xlabel = "y (km)", ylabel = "z (m)",title = "PO₄/(PO₄+kₚₒ₄)", aspect = 1)
hm_p_lim = heatmap!(ax_p_lim, yt/1e3, zt, p_lim; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_can[2, 4], hm_p_lim; flipaxis = false)
# Fe limitation
ax_f_lim = Axis(fig_can[1, 3]; xlabel = "y (km)", ylabel = "z (m)",title = "Fe/(Fe+kᵢᵣₒₙ)", aspect = 1)
hm_f_lim = heatmap!(ax_f_lim, yt/1e3, zt, f_lim; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_can[1, 4], hm_f_lim; flipaxis = false)

# All limitations
ax_all_lim = Axis(fig_can[2, 5]; xlabel = "y (km)", ylabel = "z (m)",title = "All limitations", aspect = 1)
hm_all_lim = heatmap!(ax_all_lim, yt/1e3, zt[190:200], all_lim[:,190:200]; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_can[2, 6], hm_all_lim; flipaxis = false)

display(fig_can)
=#

################################## Plot prod/remin rates ##################################


############### Calculate NCP rates ###############
# μ_surf = Field{Nothing, Center, Center}(grid) 
# μ(y,z)=  2e-5 + sinpi(y/Ly) * 8e-5  # mol P m⁻³ d⁻¹
# set!(μ_surf, μ)   
# fill_halo_regions!(μ_surf, arch)
# NCP_final = μ_surf[1,1:500,1:200] .* all_lim .* 1e3  # mmol P m⁻³ d⁻¹
# NCP_final = model.auxiliary_fields.NCP.data[1,1:500,1:200] .*1e3 .* 1days

############### Calculate Remin rates ###############
z₀ = log(0.01)*25 
rₛₑ = 0.5 # Sedimentary remineralization constant (/day)
wₛₙₖ = 10  # m/day
Δz = 20  # Vertical grid spacing in meters

# POP_r = -0.84 .* wₛₙₖ ./(z_matrix' .+ z₀) # Remineralization constant
# POP_r[:,1] .= rₛₑ # sedimentary remin in the bottom layer

# DOP_r = 2/365.25 # /day

# POP_remin_final = POP_r .* POP_final .*1e3 # mmol P m⁻³ d⁻¹
# DOP_remin_final =  DOP_r.* DOP_final .*1e3 # mmol P m⁻³ d⁻¹
# POP_remin_final = model.auxiliary_fields.Premin.data[1,1:500,1:200] .* 1e3 .* 1days
# DOP_remin_final = model.auxiliary_fields.Dremin.data[1,1:500,1:200].* 1e3 .* 1days

total_remin_final = POP_remin_final .+ DOP_remin_final

################### Plot ###################
fig_rate = Figure(size = (1300, 700))

ax_NCP = Axis(fig_rate[1, 1]; xlabel = "y (km)", ylabel = "z (m)",title = "Surface NCP (mmol m⁻³ d⁻¹)", aspect = 1)
hm_NCP = heatmap!(ax_NCP, yt/1e3, zt[190:200], NCP_final[:,190:200]; colorrange = (0,5e-3),colormap = :viridis) 
Colorbar(fig_rate[1, 2], hm_NCP; flipaxis = false)

ax_remin = Axis(fig_rate[2, 1]; xlabel = "y (km)", ylabel = "z (m)",title = "POP remin (mmol m⁻³ d⁻¹)", aspect = 1)
hm_remin = heatmap!(ax_remin, yt/1e3, zt, POP_remin_final; colorrange = (0,4e-4),colormap = :viridis) 
Colorbar(fig_rate[2, 2], hm_remin; flipaxis = false)

total_NCP = sum(NCP_final; dims = 2).*Δz # mmol m⁻² d⁻¹
ax_total_NCP = Axis(fig_rate[1, 4]; xlabel = "y (km)", ylabel = "Rate (mmol m⁻² d⁻¹)", title = "Integrated NCP")
ylims!(ax_total_NCP, 0,0.5)
lines!(ax_total_NCP, yt/1e3, vec(total_NCP),linewidth = 2)

total_R = sum(total_remin_final; dims = 2).*Δz 
total_PR = sum(POP_remin_final; dims = 2).*Δz 
total_DR = sum(DOP_remin_final; dims = 2).*Δz 
ax_total_R = Axis(fig_rate[2, 4]; xlabel = "y (km)", ylabel = "Rate (mmol m⁻² d⁻¹)", title = "Integrated Remin")
ylims!(ax_total_R, 0, 0.5)
lines!(ax_total_R, yt/1e3, vec(total_PR),linewidth = 2, label = "POP")
lines!(ax_total_R, yt/1e3, vec(total_DR),linewidth = 2, label = "DOP")
lines!(ax_total_R, yt/1e3, vec(total_R),linewidth = 2, label = "Total")
axislegend(ax_total_R, position = :rt)

# ax_reminD = Axis(fig_rate[3, 1]; xlabel = "y (km)", ylabel = "z (m)",title = "DOP remin (mmol m⁻³ d⁻¹)", aspect = 1)
# hm_reminD = heatmap!(ax_reminD, yt/1e3, zt, DOP_remin_final; colorrange = (0,2e-3),colormap = :viridis) 
# Colorbar(fig_rate[3, 2], hm_reminD; flipaxis = false)

avg_NCP = mean(NCP_final; dims = 1)
ax_avg_NCP = Axis(fig_rate[1, 3]; xlabel = "Rate (mmol m⁻³ d⁻¹)", ylabel = "z (m)", title = "Ave. NCP")
xlims!(ax_avg_NCP, -5e-4, 5e-3)
lines!(ax_avg_NCP, vec(avg_NCP),zt, linewidth = 2)

avg_PR = mean(POP_remin_final; dims = 1)
avg_DR = mean(DOP_remin_final; dims = 1)
ax_avg_R = Axis(fig_rate[2, 3]; xlabel = "Rate (mmol m⁻³ d⁻¹)", ylabel = "z (m)", title = "Ave. Remin")
xlims!(ax_avg_R, -1e-6, 1.1e-3)
lines!(ax_avg_R, vec(avg_PR),zt, linewidth = 2, label = "POP")
lines!(ax_avg_R, vec(avg_DR),zt, linewidth = 2, label = "DOP")
axislegend(ax_avg_R, position = :rb)

# Calculate R/NCP as a function of latitude
RPratio = total_R ./ total_NCP # mmol m⁻² d⁻¹

ax_RP = Axis(fig_rate[1, 5]; xlabel = "y (km)", ylabel = "ratio", title = "Remin / NCP")
ylims!(ax_RP, 0.9, 1.3)
# Remove boundaries: 0 PAR -> 0 NCP
lines!(ax_RP, yt[5:495]/1e3, vec(RPratio[5:495]), linewidth = 2)

# Calculate e-ratio 
eratio_100 = wₛₙₖ .* POP_final[:,196] .* 1e3 ./ total_NCP # mmol m⁻² d⁻¹

ax_eR = Axis(fig_rate[2, 5]; xlabel = "y (km)", ylabel = "e-ratio", title = "e-ratio (F₁₀₀ₘ/NCP)")
ylims!(ax_eR, 0.18, 0.22)
# Remove boundaries: 0 PAR -> 0 NCP
lines!(ax_eR, yt[5:495]/1e3, vec(eratio_100[5:495]), linewidth = 2)

display(fig_rate)