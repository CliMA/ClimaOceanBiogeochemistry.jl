# Plotting codes: visualize 1-D column model results
# NPD vs NPBD models

# Load data
using Oceananigans
using GLMakie
using JLD2

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

Nn_last = interior(Nt[end], 1, 1, :)
Pn_last = interior(Pt[end], 1, 1, :)
Bn_last = interior(Bt[end], 1, 1, :)
Dn_last = interior(Dt[end], 1, 1, :)

N1n_last = interior(N1t[end], 1, 1, :)
P1n_last = interior(P1t[end], 1, 1, :)
D1n_last = interior(D1t[end], 1, 1, :)

# Start plotting
fig = Figure(;size=(1200, 800))

# NPBD model
axN = Axis(fig[1, 1], ylabel="z (m)", xlabel="[N] (mmol m⁻³)")
xlims!(axN, -0.1, 50)
axB = Axis(fig[1, 2], xlabel="[P] or [B] (mmol m⁻³)")
xlims!(axB, -0.05, 0.75)
axD = Axis(fig[1, 3], xlabel="[D] (mmol m⁻³)")
xlims!(axD, -0.1, 1.5)

lines!(axN, Nn_last, z)
lines!(axB, Pn_last, z, label="P")
lines!(axB, Bn_last, z, label="B")
lines!(axD, Dn_last, z)

axislegend(axB, position = :rb)

# NPD model
axN1 = Axis(fig[2, 1], ylabel="z (m)", xlabel="[N] (mmol m⁻³)")
axB1 = Axis(fig[2, 2], xlabel="[P] (mmol m⁻³)")
axD1 = Axis(fig[2, 3], xlabel="[D] (mmol m⁻³)")

lines!(axN1, N1n_last, z)
lines!(axB1, P1n_last, z)
lines!(axD1, D1n_last, z)

display(fig)
# Save the last frame as a figure
save("test.png", fig)
