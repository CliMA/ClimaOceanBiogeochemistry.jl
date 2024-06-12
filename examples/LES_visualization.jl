using Statistics
using GLMakie
using Oceananigans
using Oceananigans.Units
using JLD2

filename = "LES_test4"
data = load(filename * ".jld2")

# Output
Nt = FieldTimeSeries(filename * ".jld2", "N")
Nt_end = Nt[:,:,:,end]
Pt = FieldTimeSeries(filename * ".jld2", "P")
Pt_end = Pt[:,:,:,end]
Bt = FieldTimeSeries(filename * ".jld2", "B")
Bt_end = Bt[:,:,:,end]
Dt = FieldTimeSeries(filename * ".jld2", "D")
Dt_end = Dt[:,:,:,end]

Nn = interior(Nt_end, :, 1, :)
Pn = interior(Pt_end, :, 1, :)
Bn = interior(Bt_end, :, 1, :)
Dn = interior(Dt_end, :, 1, :)

fig1 = Figure()

axN = Axis(fig1[1, 1], title="Nutrients")
axP = Axis(fig1[1, 2], title="Phytoplankton")
axB = Axis(fig1[1, 3], title="Bacteria")
axD = Axis(fig1[1, 4], title="Detritus")

heatmap!(axN, Nn)
heatmap!(axP, Pn)
heatmap!(axB, Bn)
heatmap!(axD, Dn)

N_col = vec(mean(Nn, dims=1))
P_col = vec(mean(Pn, dims=1))
B_col = vec(mean(Bn, dims=1))
D_col = vec(mean(Dn, dims=1))

axNcol = Axis(fig1[2, 1], title="Nutrients")
axPcol = Axis(fig1[2, 2], title="Phytoplankton")
axBcol = Axis(fig1[2, 3], title="Bacteria")
axDcol = Axis(fig1[2, 4], title="Detritus")

z = vec(reshape(collect(1.0:64.0), 1, 64))
lines!(axNcol, N_col, z)
lines!(axPcol, P_col, z)
lines!(axBcol, B_col, z)
lines!(axDcol, D_col, z)

display(fig1)
save("NPZDB_4.png", fig1)
