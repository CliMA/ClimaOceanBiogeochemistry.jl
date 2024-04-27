# Plotting codes: visualize 1-D column model results
# NPD vs NPBD models
# Rates: 

# Load data
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: znode,Center, AbstractTopology, Flat, Bounded
const c = Center()
using GLMakie

# filename1 = "plot/NPBD_3.jld2"
# filename2 = "plot/NPD_2.jld2"
filename1 = "plot/pert_NPBD_10y.jld2"
filename2 = "plot/pert_NPD_10y.jld2"

Nt = FieldTimeSeries(filename1, "N")
Pt = FieldTimeSeries(filename1, "P")
Bt = FieldTimeSeries(filename1, "B")
Dt = FieldTimeSeries(filename1, "D")

N1t = FieldTimeSeries(filename2, "N")
P1t = FieldTimeSeries(filename2, "P")
D1t = FieldTimeSeries(filename2, "D")

t = Nt.times
nt = length(t)
z = znodes(Nt)

# Calculate rates: primary production and respiration 

# pp = mu_P*P (phytoplankton production)
μᵖ = 1 # original unit s⁻¹ but I plot d⁻¹
kᴺ = 0.1
kᴵ = 10
λ = 25
I = 700 * exp.(z / λ)
@inline R_production(N, P) = μᵖ .* min.(N ./ (N .+ kᴺ) , I ./ (I .+ kᴵ)) .* P

# respiration = (1/y - 1)*mu_B*B (production of inorganic nutrient again by heterotrophs).
# Nonlinear (involing B)
μᵇ = 0.9
kᴰ = 0.05
y = 0.2
@inline R_respiration(D, B) = (1-y) .* μᵇ .* D ./ (D .+ kᴰ) .* B 
# Linear (no B)
r = 0.1
@inline R_remineralization(D) = r .* D

# each variable vs. time
N_time = zeros(length(t),length(z))
P_time = zeros(length(t),length(z))
B_time = zeros(length(t),length(z))
D_time = zeros(length(t),length(z))
Rate_PP = zeros(length(t),length(z))
Rate_R = zeros(length(t),length(z))

N1_time = zeros(length(t),length(z))
P1_time = zeros(length(t),length(z))
D1_time = zeros(length(t),length(z))
Rate_PP1 = zeros(length(t),length(z))
Rate_R1 = zeros(length(t),length(z))

for i in 1:nt
    N_time[i,:] = interior(Nt[i], 1, 1, :)
    P_time[i,:] = interior(Pt[i], 1, 1, :)
    B_time[i,:] = interior(Bt[i], 1, 1, :)
    D_time[i,:] = interior(Dt[i], 1, 1, :)

    Rate_PP[i,:] = R_production(N_time[i,:], P_time[i,:])
    Rate_R[i,:] = R_respiration(D_time[i,:], B_time[i,:])

    N1_time[i,:] = interior(N1t[i], 1, 1, :)
    P1_time[i,:] = interior(P1t[i], 1, 1, :)
    D1_time[i,:] = interior(D1t[i], 1, 1, :)

    Rate_PP1[i,:] = R_production(N1_time[i,:], P1_time[i,:])
    Rate_R1[i,:] = R_remineralization(D1_time[i,:])
end

# The end point
Rate_PP_last = Rate_PP[end,:]
Rate_R_last = Rate_R[end,:]

Rate_PP1_last = Rate_PP1[end,:]
Rate_R1_last = Rate_R1[end,:]

fig = Figure(;size=(1500, 1500))
# Row 1: NPBD
axPP = Axis(fig[1, 1])
hmPP = heatmap!(axPP, Rate_PP, colormap=:plasma, colorrange = (0, 0.8))
axPP.xlabel = "Time"
axPP.ylabel = "Depth"
axPP.title = " Primary production rate (mmol m⁻³ d⁻¹)"
Colorbar(fig[1, 2], hmPP)
 
axR = Axis(fig[1, 3])
hmR = heatmap!(axR, Rate_R, colormap=:plasma, colorrange = (0, 0.1))
axR.xlabel = "Time"
axR.ylabel = "Depth"
axR.title = " Respiration rate (mmol m⁻³ d⁻¹)"
Colorbar(fig[1, 4], hmR)

# Define custom colormap
function custom_colormap(x)
    if x <= 1
        return HSL(240, 1, 0.5*(1+x))  # blue to white transition (Hue 240 is blue)
    else
        return HSL(0, 1, 1-0.5*(x-1)/9)  # white to red transition (Hue 0 is red)
    end
end

cmap = [custom_colormap(x) for x in range(0, stop=10, length=256)]

axRatio = Axis(fig[1, 5])
hmRatio = heatmap!(axRatio, Rate_PP./Rate_R, colormap=cmap, colorrange = (0, 10))
axRatio.xlabel = "Time"
axRatio.ylabel = "Depth"
axRatio.title = " Production : Respiration"
Colorbar(fig[1, 6], hmRatio)

# Row 2: NPD
axPP1 = Axis(fig[2, 1])
hmPP1 = heatmap!(axPP1, Rate_PP1, colormap=:plasma, colorrange = (0, 0.8))
axPP1.xlabel = "Time"
axPP1.ylabel = "Depth"
axPP1.title = " Primary production rate (mmol m⁻³ d⁻¹)"
Colorbar(fig[2, 2], hmPP1)
 
axR1 = Axis(fig[2, 3])
hmR1 = heatmap!(axR1, Rate_R1, colormap=:plasma, colorrange = (0, 0.1))
axR1.xlabel = "Time"
axR1.ylabel = "Depth"
axR1.title = " Respiration rate (mmol m⁻³ d⁻¹)"
Colorbar(fig[2, 4], hmR1)

axRatio1 = Axis(fig[2, 5])
hmRatio1 = heatmap!(axRatio1, Rate_PP1./Rate_R1, colormap=cmap, colorrange = (0, 10))
axRatio1.xlabel = "Time"
axRatio1.ylabel = "Depth"
axRatio1.title = " Production : Respiration"
Colorbar(fig[2, 6], hmRatio1)

# Row 3: vertical profiles
axPP2 = Axis(fig[3, 1], ylabel="z (m)", xlabel="Primary production rate (mmol m⁻³ d⁻¹)")
axR2 = Axis(fig[3, 3], xlabel="Respiration rate (mmol m⁻³ d⁻¹)")
axRatio2 = Axis(fig[3, 5], xlabel=" Ratio of production to respiration")

lines!(axPP2, Rate_PP_last, z, label="NPBD")
lines!(axPP2, Rate_PP1_last, z, label="NPD")
axislegend(axPP2, position = :rb)

lines!(axR2, Rate_R_last, z, label="NPBD")
lines!(axR2, Rate_R1_last, z, label="NPD")
axislegend(axR2, position = :rb)

lines!(axRatio2, Rate_PP_last./Rate_R_last, z, label="NPBD")
lines!(axRatio2, Rate_PP1_last./Rate_R1_last, z, label="NPD")
axislegend(axRatio2, position = :rb)

display(fig)
save("test.png", fig)
