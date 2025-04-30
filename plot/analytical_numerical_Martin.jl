using GLMakie
using Oceananigans
using Statistics

z₀ = log(0.01)*25 

# Create z and t ranges
z_range = -1000:20:0   
t_range = 1:1:730 


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
    temp1 = 1e4*interior(POP_timeseries1[i], 1, :, :)
    POP_1[i, :] = mean(temp1,dims=1)
    temp2 = 1e4*interior(POP_timeseries2[i], 1, :, :)
    POP_2[i, :] = mean(temp2,dims=1)
 end
 for i in 366:730
    temp1 = 1e4*interior(POP_timeseries1[i-365], 1, :, :)
    POP_1[i, :] = mean(temp1,dims=1)
    temp2 = 1e4*interior(POP_timeseries2[i-365], 1, :, :)
    POP_2[i, :] = mean(temp2,dims=1)
 end

POP_ref1 = POP_1[:, 195]
POP_ref2 = POP_2[:, 195]

using Interpolations
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


# Define the function D(z, t)
# z and t are positional argument that must be provided
# the other are default numbers that can be changed

# function D1(z, t; D0=1.0, A = 1.0, N=365, w=-10, z0=z₀, b=0.84)
#     return (D0+A*sinpi(2 / N * (t - (z-z0*2)/w))) * ((z + z0) / (2* z0))^(-b) 
# end
# function D2(z, t; D0=1.0, A = 1.0, N=365/12, w=-10, z0=z₀, b=0.84)
#     return (D0+A*sinpi(2 / N * (t - (z-z0/2)/w))) * ((z + z0) / (2* z0))^(-b) 
# end


matrix1 = compute_D(collect(t_range), POP_ref1, collect(zw), z₀, 0.84,-10.0)
matrix2 = compute_D(collect(t_range), POP_ref2, collect(zw), z₀, 0.84,-10.0)

z₀ = log(0.01)*25 
Martin_factor = ((2*z₀) ./ (zw .+ z₀)).^0.84 # zw[190] .+ 

# POP flux in Martin curve (mmol m-2 d-1)
avg_MartinFlux1 = zeros(730, 200) 
avg_MartinFlux2 = zeros(730, 200) 
for i in 1:730
    temp_Martin1 = repeat(POP_ref1[i,:], 1, 200)
    temp_Martin2 = repeat(POP_ref2[i,:], 1, 200)
    for k in 1:200
        temp_Martin1[:,k] .= temp_Martin1[:,k] .* Martin_factor[k]
        temp_Martin2[:,k] .= temp_Martin2[:,k] .* Martin_factor[k]
    end
    # Calculate domain-average flux 
    avg_MartinFlux1[i,:] = temp_Martin1
    avg_MartinFlux2[i,:] = temp_Martin2
end
#=
D0_1 = mean(POP_ref1)
D0_2 = mean(POP_ref2)
A_1 = maximum(POP_ref1)-mean(POP_ref1)#(maximum(POP_ref1)-minimum(POP_ref1))/2
A_2 = maximum(POP_ref2)-mean(POP_ref2)#(maximum(POP_ref2)-minimum(POP_ref2))/2
N_1 = 365
N_2 = 365/12

# Evaluate D(z, t) over the grid
D1_values = [D1(z, t; D0=D0_1, A = A_1, N=N_1) for z in zw, t in t_range]
D2_values = [D2(z, t; D0=D0_2, A = A_2, N=N_2) for z in zw, t in t_range]
matrix1 = [D1_values[i, j][1] for j in 1:365, i in 1:200]
matrix2 = [D2_values[i, j][1] for j in 1:365, i in 1:200]
=#

# Plot
#=
fig = Figure()

ax_t_eR = Axis(fig[1, 1];  xlabel = "t (day)", title = "yearly",
            yticklabelcolor=:brown3, ylabelcolor=:brown3)
lines!(ax_t_eR, 1:1:365, POP_ref1, color=:brown3, linewidth = 2, label = "N")
lines!(ax_t_eR, 1:1:365, matrix1'[:, 195], color=:black, linewidth = 2, label = "A")
axislegend(ax_t_eR, position = :rt)

ax_t_2 = Axis(fig[1, 2];  xlabel = "t (day)", title = "monthly",
            yticklabelcolor=:brown3, ylabelcolor=:brown3)
lines!(ax_t_2, 1:1:365, POP_ref2, color=:brown3, linewidth = 2, label = "N")
lines!(ax_t_2, 1:1:365, matrix2'[:, 195], color=:black, linewidth = 2, label = "A")
axislegend(ax_t_2, position = :lt)

display(fig)
=#

# convert P to C
P2Cratio = 117

xticks = 365:100:730  # or whatever spacing you want
xtick_labels = string.(xticks .- 365)

fig = Figure()

##### Analytical #####
ax1 = Axis(fig[1, 1], xlabel="t (days)", ylabel="Depth (m)", title="Analytical POC (mmol m⁻³ d⁻¹)")
hm_1 = heatmap!(ax1, t_range, zw, P2Cratio .* matrix1', colorrange = (0, 8),colormap=:thermal, interpolate = true)
Colorbar(fig[1, 2], hm_1; flipaxis = false)
ylims!(ax1, -1000,-100)
xlims!(ax1, 366, 730)
ax1.xticks = (xticks, xtick_labels)

ax2 = Axis(fig[2, 1], xlabel="t (days)", ylabel="Depth (m)", title="Analytical POC (mmol m⁻³ d⁻¹)")
hm_2 = heatmap!(ax2, t_range, zw, P2Cratio .* matrix2', colorrange = (0, 8),colormap=:thermal, interpolate = true)
Colorbar(fig[2, 2], hm_2; flipaxis = false)
ylims!(ax2, -1000,-100)
xlims!(ax2, 366, 730)
ax2.xticks = (xticks, xtick_labels)

##### Numerical #####
ax_True1 = Axis(fig[1, 3]; xlabel = "t (days)", ylabel = "Depth (m)", title = "Numerical POC (mmol m⁻³ d⁻¹)", aspect = 1)
hm_True1 = heatmap!(ax_True1, t_range, zw, P2Cratio .* POP_1; colorrange = (0, 8),colormap = :thermal, interpolate = true) 
Colorbar(fig[1, 4], hm_True1; flipaxis = false)
ylims!(ax_True1, -1000, -100)
xlims!(ax_True1, 366, 730)
ax_True1.xticks = (xticks, xtick_labels)

ax_True2 = Axis(fig[2, 3]; xlabel = "t (days)", ylabel = "Depth (m)", title = "Numerical POC (mmol m⁻³ d⁻¹)", aspect = 1)
hm_True2 = heatmap!(ax_True2, t_range, zw, P2Cratio .* POP_2; colorrange = (0, 8),colormap = :thermal, interpolate = true) 
Colorbar(fig[2, 4], hm_True2; flipaxis = false)
ylims!(ax_True2, -1000, -100)
xlims!(ax_True2, 366, 730)
ax_True2.xticks = (xticks, xtick_labels)

############# Martin heatmap #############
ax_Martin1 = Axis(fig[1, 5]; xlabel = "t (days)", ylabel = "Depth (m)", title = "Martin POC (mmol m⁻³ d⁻¹)", aspect = 1)
hm_Martin1 = heatmap!(ax_Martin1, t_range, zw, P2Cratio .* avg_MartinFlux1; colorrange = (0, 8),colormap = :thermal, interpolate = true) 
Colorbar(fig[1, 6], hm_Martin1; flipaxis = false)
ylims!(ax_Martin1, -1000, -100)
xlims!(ax_Martin1, 366, 730)
ax_Martin1.xticks = (xticks, xtick_labels)

ax_Martin2 = Axis(fig[2, 5]; xlabel = "t (days)", ylabel = "Depth (m)", title = "Martin POC (mmol m⁻³ d⁻¹)", aspect = 1)
hm_Martin2 = heatmap!(ax_Martin2, t_range, zw, P2Cratio .* avg_MartinFlux2; colorrange = (0, 8),colormap = :thermal, interpolate = true) 
Colorbar(fig[2, 6], hm_Martin2; flipaxis = false)
ylims!(ax_Martin2, -1000, -100)
xlims!(ax_Martin2, 366, 730)
ax_Martin2.xticks = (xticks, xtick_labels)

display(fig)


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


#
#=

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