using GLMakie
using Oceananigans
using Statistics
using Interpolations
using XLSX

xlsx_data = XLSX.readxlsx("bats_flux.xlsx")
OC = xlsx_data["reconstruct"]

dayofyear = OC[2:21, 3]  # after 2020
OCflux_150 = OC[2:21, 4]
OCflux_300 = OC[2:21, 5]

z₀ = -110.0 #log(0.01)*25 

# Create z and t ranges
z_range = collect(-1000:10:0)   
t_range = vec(dayofyear)

# The boundary condition at z0 is defined as D0, which is a vector of times
function compute_D(t::Vector, D0::Vector, z::Vector, z0::Float64, b::Float64, w::Float64)
    t_grid = reshape(t, 1, :)   # (1, Nt)
    z_grid = reshape(z, :, 1)   # (Nz, 1)
    t_shifted = clamp.(t_grid .- (z_grid .- z0) ./ w, minimum(t), maximum(t))
    interp_D0 = linear_interpolation(t, D0)
    F_shifted = interp_D0.(t_shifted)
    Martin_term = ((z .+ z0) ./ (2*z0)).^(-b)
    return F_shifted .* Martin_term
end

# P_ref = 1 .+ 0.5 .*sinpi.(2 .* t_range ./ 365*2)
# P_ref = zeros(length(t_range)) 
# P_ref[10:30] .= 1.0
# P_ref[10:39] .= sinpi.((t_range[10:39]) ./ 30)
P = compute_D(t_range, vec(OCflux_150), collect(z_range), -150.0, 0.84,-10.0)

# Average with time
P_avg = mean(P; dims=2)

# Plotting
# fig = Figure()

ax0 = Axis(fig[2, 1]; xlabel = "t (days)", ylabel = "Concentration", title = "Particle export out of z₀", aspect = 1)
lines!(ax0, t_range, vec(OCflux_150), color=:black, linewidth = 3)
# ylims!(ax0, 0, 1.1)
# xlims!(ax0, 0, 150)

ax1 = Axis(fig[2, 2]; xlabel = "t (days)", ylabel = "Depth (m)", title = "Propagation and attenuation of particle pulse", aspect = 1)
hm1 = heatmap!(ax1, t_range, z_range, P'; colorrange = (0, 80),colormap = :bluegreenyellow, interpolate = true) 
Colorbar(fig[2, 3], hm1; flipaxis = false)
ylims!(ax1, -1000, -150)
# xlims!(ax1, 0, 150)


z_ind = [81, 41, 1]     # pick a few depths: 200, 800, 1000
t_ind = [4, 12, 19]  # pick a few time snapshots: 210, 460, 685

# --- Left Panel: Vertical profiles at selected times ---
#
ax11 = Axis(fig[2, 4], xlabel="Particle", ylabel="Depth (m)", title="Vertical particle profiles at different times")

lines!(ax11, P[:, t_ind[1]], z_range; label = "t = $(t_range[t_ind[1]]) days", linewidth = 3)
lines!(ax11, P[:, t_ind[2]], z_range; label = "t = $(t_range[t_ind[2]]) days", linewidth = 3)
lines!(ax11, P[:, t_ind[3]], z_range; label = "t = $(t_range[t_ind[3]]) days", linewidth = 3)
lines!(ax11, vec(P_avg), z_range; label = "2-year average", color=:black, linestyle=:dash, linewidth = 3)
axislegend(ax11, position=:rb)
ylims!(ax11, -1000, 0)
# xlims!(ax11, 0, 1)

# --- Right Panel: Time series at selected depths ---
ax2 = Axis(fig[2, 5], xlabel="Time (days)", ylabel="Concentration", title="Particle time series at fixed depths")
lines!(ax2, t_range, P[z_ind[1], :]; label="z = $(z_range[z_ind[1]]) m", linewidth = 2)
lines!(ax2, t_range, P[z_ind[2], :]; label="z = $(z_range[z_ind[2]]) m", linewidth = 2)
lines!(ax2, t_range, P[z_ind[3], :]; label="z = $(z_range[z_ind[3]]) m", linewidth = 2)
axislegend(ax2, position=:lt)
# xlims!(ax2, 0, 150)
#
display(fig)
#
#=
##################### Load model #####################
filepath1 = "./P5_51y.jld2"
filepath2 = "./P6_51y.jld2"

POP_timeseries1 = FieldTimeSeries(filepath1, "POP")
POP_timeseries2 = FieldTimeSeries(filepath2, "POP")
times = POP_timeseries1.times
xw, yw, zw = nodes(POP_timeseries1)

POP_1 = zeros(730, 200) 
POP_2 = zeros(730, 200) 
for i in 1:365
    temp1 = 1e3*interior(POP_timeseries1[i], 1, :, :)
    POP_1[i, :] = mean(temp1,dims=1)
    temp2 = 1e3*interior(POP_timeseries2[i], 1, :, :)
    POP_2[i, :] = mean(temp2,dims=1)
 end
 for i in 366:730
    temp1 = 1e3*interior(POP_timeseries1[i-365], 1, :, :)
    POP_1[i, :] = mean(temp1,dims=1)
    temp2 = 1e3*interior(POP_timeseries2[i-365], 1, :, :)
    POP_2[i, :] = mean(temp2,dims=1)
 end

POP_ref1 = POP_1[:, 195]
POP_ref2 = POP_2[:, 195]

matrix1 = compute_D(collect(t_range), POP_ref1, collect(zw), z₀, 0.84,-10.0)
matrix2 = compute_D(collect(t_range), POP_ref2, collect(zw), z₀, 0.84,-10.0)

fig = Figure()

ax_True1 = Axis(fig[1, 1]; xlabel = "t (days)", ylabel = "Depth (m)", title = "Numerical POP (mmol m⁻³ d⁻¹)", aspect = 1)
hm_True1 = heatmap!(ax_True1, t_range, zw, POP_1; colorrange = (0, 0.009),colormap = :thermal, interpolate = true) 
Colorbar(fig[1, 2], hm_True1; flipaxis = false)
ylims!(ax_True1, -1000, -100)
xlims!(ax_True1, 366, 730)

ax_True2 = Axis(fig[2, 1]; xlabel = "t (days)", ylabel = "Depth (m)", title = "Numerical POP (mmol m⁻³ d⁻¹)", aspect = 1)
hm_True2 = heatmap!(ax_True2, t_range, zw, POP_2; colorrange = (0, 0.009),colormap = :thermal, interpolate = true) 
Colorbar(fig[2, 2], hm_True2; flipaxis = false)
ylims!(ax_True2, -1000, -100)
xlims!(ax_True2, 366, 730)

ax1 = Axis(fig[1, 3], xlabel="t (days)", ylabel="Depth (m)", title="Analytical D(z, t)")
hm_1 = heatmap!(ax1, t_range, zw, matrix1', colorrange = (0, 0.009),colormap=:thermal, interpolate = true)
Colorbar(fig[1, 4], hm_1; flipaxis = false)
ylims!(ax1, -1000,-100)
xlims!(ax1, 366, 730)

ax2 = Axis(fig[2, 3], xlabel="t (days)", ylabel="Depth (m)", title="Analytical D(z, t)")
hm_2 = heatmap!(ax2, t_range, zw, matrix2', colorrange = (0, 0.009),colormap=:thermal, interpolate = true)
Colorbar(fig[2, 4], hm_2; flipaxis = false)
ylims!(ax2, -1000,-100)
xlims!(ax2, 366, 730)

# ax1diff = Axis(fig[1, 5], xlabel="t (days)", ylabel="Depth (m)", title="Analytical - Numerical")
# hm_1diff = heatmap!(ax1diff, t_range, zw, (matrix1'.-POP_1), colorrange = (-0.002, 0.002), colormap=:balance, interpolate = true)
# Colorbar(fig[1, 6], hm_1diff; flipaxis = false)
# ylims!(ax1diff, -1000,-100)
# xlims!(ax1diff, 366, 730)

# ax2diff = Axis(fig[2, 5], xlabel="t (days)", ylabel="Depth (m)", title="Analytical - Numerical")
# hm_2diff = heatmap!(ax2diff, t_range, zw, (matrix2'.-POP_2), colorrange = (-0.002, 0.002), colormap=:balance, interpolate = true)
# Colorbar(fig[2, 6], hm_2diff; flipaxis = false)
# ylims!(ax2diff, -1000,-100)
# xlims!(ax2diff, 366, 730)

display(fig)

# Teff
z_factor = ((-1010+z₀)/(-210+z₀))^-0.84
phase1 = (-1010-z₀/2)/(-10)
phase2 = (-210-z₀/2)/(-10)

T_eff1 = (D0_1.+A_1.*sinpi.(2/N_1.*(t_range.-phase1)))./(D0_1.+A_1.*sinpi.(2/N_1.*(t_range.-phase2))).*z_factor
T_eff2 = (D0_2.+A_2.*sinpi.(2/N_2.*(t_range.-phase1)))./(D0_2.+A_2.*sinpi.(2/N_2.*(t_range.-phase2))).*z_factor

POP_Teff1 = POP_1[:, 150] ./ POP_1[:, 190]
POP_Teff2 = POP_2[:, 150] ./ POP_2[:, 190]

fig = Figure()

ax_tT = Axis(fig[1, 1];  xlabel = "t (day)", title = "Yearly: F₁₀₀₀/F₂₀₀",
            yticklabelcolor=:brown3, ylabelcolor=:brown3)
lines!(ax_tT, 1:1:365, T_eff1, color=:brown3, linewidth = 2, label = "Analytical")
lines!(ax_tT, 1:1:365, POP_Teff1, color=:black, linewidth = 2, label = "Numerical")
ylims!(ax_tT, 0.1,0.85)
axislegend(ax_tT, position = :lt)

ax_tT2 = Axis(fig[1, 2];  xlabel = "t (day)", title = "Monthly: F₁₀₀₀/F₂₀₀",
            yticklabelcolor=:brown3, ylabelcolor=:brown3)
lines!(ax_tT2, 1:1:365, T_eff2, color=:brown3, linewidth = 2, label = "Analytical")
lines!(ax_tT2, 1:1:365, POP_Teff2, color=:black, linewidth = 2, label = "Numerical")
ylims!(ax_tT2, 0.1,0.85)
axislegend(ax_tT2, position = :lt)

display(fig)
=#