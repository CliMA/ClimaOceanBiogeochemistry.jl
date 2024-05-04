# Plotting codes: visualize 1-D column model results
# NPD vs NPBD models
# Rates: 

# Load data
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: znode,Center, AbstractTopology, Flat, Bounded
const c = Center()
using GLMakie

filename1 = "plot/NPBD_3.jld2"
# filename2 = "plot/NPD_2.jld2"
# filename1 = "plot/pert_NPBD_10y.jld2"

Bt = FieldTimeSeries(filename1, "B")
Dt = FieldTimeSeries(filename1, "D")

t = Bt.times
nt = length(t)
z = znodes(Bt)

# Calculate ‘equivalent R’ in NPBD: = Vmax/(D+kD)*B
μᵇ = 0.9
kᴰ = 0.05
@inline Req(D, B) = μᵇ .* B ./ (D .+ kᴰ) 

# each variable vs. time
B_time = zeros(length(t),length(z))
D_time = zeros(length(t),length(z))
Req_time = zeros(length(t),length(z))

for i in 1:nt
    B_time[i,:] = interior(Bt[i], 1, 1, :)
    D_time[i,:] = interior(Dt[i], 1, 1, :)

    Req_time[i,:] = Req(D_time[i,:], B_time[i,:])
end
Req_last = Req_time[end,:]

fig = Figure(;size=(500, 500))

ax = Axis(fig[1, 1], ylabel="z (m)", xlabel="Equivalent r")

r = zeros(length(z)).+0.1
lines!(ax, Req_last, z, label="req in NPBD")
lines!(ax, r, z, label="r in NPD")
axislegend(ax, position = :rb)

display(fig)
save("test.png", fig)
