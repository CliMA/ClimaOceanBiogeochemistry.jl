using Statistics
using GLMakie
using Oceananigans.Units
using JLD2

# Output
data = JLD2.load("LES_test_xz.jld2")
#=
variable1 = data["variable1"]
variable2 = data["variable2"]
variable3 = data["variable3"]
=#

# data 
T = model.tracers.T
N = model.tracers.N
P = model.tracers.P
Z = model.tracers.Z
B = model.tracers.B
D = model.tracers.D
w = model.velocities.w

Tn = interior(T, :, 1, :)
wn = interior(w, :, 1, :)
Nn = interior(N, :, 1, :)
Pn = interior(P, :, 1, :)
#Zn = interior(Z, :, 1, :)
Bn = interior(B, :, 1, :)
Dn = interior(D, :, 1, :)

fig1 = Figure()

axT = Axis(fig1[1, 1], title="Temperature")
axw = Axis(fig1[1, 2], title="Vertical velocity")
axN = Axis(fig1[1, 3], title="Nutrients")
axP = Axis(fig1[2, 1], title="Phytoplankton")
#axZ = Axis(fig[2, 2], title="Zooplankton")
axB = Axis(fig1[2, 2], title="Bacteria")
axD = Axis(fig1[2, 3], title="Detritus")

heatmap!(axT, Tn)
heatmap!(axw, wn)
heatmap!(axN, Nn)
heatmap!(axP, Pn)
#heatmap!(axZ, Zn)
heatmap!(axB, Bn)
heatmap!(axD, Dn)

display(fig1)
save("NPZDB_hm5.png", fig1)

fig2 = Figure()
N_col = vec(mean(Nn, dims=1))
P_col = vec(mean(Pn, dims=1))
B_col = vec(mean(Bn, dims=1))
D_col = vec(mean(Dn, dims=1))

axNcol = Axis(fig2[1, 1], title="Nutrients")
axPcol = Axis(fig2[1, 2], title="Phytoplankton")
axBcol = Axis(fig2[2, 1], title="Bacteria")
axDcol = Axis(fig2[2, 2], title="Detritus")

z = vec(reshape(collect(1.0:64.0), 1, 64))
lines!(axNcol, N_col, z)
lines!(axPcol, P_col, z)
lines!(axBcol, B_col, z)
lines!(axDcol, D_col, z)

display(fig2)
save("NPZDB_line5.png", fig2)