# Plotting codes: visualize 1-D column model results
# NPD vs NPBD models

# Load data
using Oceananigans
using GLMakie
using JLD2

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

# Figure : sum of each variable vs. time
N_time = zeros(1:nt)
P_time = zeros(1:nt)
B_time = zeros(1:nt)
D_time = zeros(1:nt)

N1_time = zeros(1:nt)
P1_time = zeros(1:nt)
D1_time = zeros(1:nt)

for times = 1:nt
    N_time[times] = sum(Nt[:,:,:,times])
    P_time[times] = sum(Pt[:,:,:,times])
    B_time[times] = sum(Bt[:,:,:,times])
    D_time[times] = sum(Dt[:,:,:,times])
    N1_time[times] = sum(N1t[:,:,:,times])
    P1_time[times] = sum(P1t[:,:,:,times])
    D1_time[times] = sum(D1t[:,:,:,times])
end

TimeVar = Figure()
ax1 = Axis(TimeVar[1,1], title="Variable over time",ylabel="Variable (mmol m⁻³)", xlabel="Time (days)")
lines!(ax1, 1:nt, N_time, label="N")  
lines!(ax1, 1:nt, P_time, label="P")  
lines!(ax1, 1:nt, B_time, label="B")  
lines!(ax1, 1:nt, D_time, label="D") 
axislegend()
 
ax2 = Axis(TimeVar[2,1], title="Variable over time",ylabel="Variable (mmol m⁻³)", xlabel="Time (days)")
lines!(ax2, 1:nt, N1_time, label="N")  
lines!(ax2, 1:nt, P1_time, label="P")   
lines!(ax2, 1:nt, D1_time, label="D") 
axislegend()

save("test.png", TimeVar)
