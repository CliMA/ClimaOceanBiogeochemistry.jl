# Plotting codes: visualize 1-D column model results
# NPD vs NPBD models

# Load data
using Oceananigans
using GLMakie

filename1 = "plot/NPBDF_2.jld2"
filename3 = "plot/NPDF_2.jld2"

Nt = FieldTimeSeries(filename1, "N")
Pt = FieldTimeSeries(filename1, "P")
Bt = FieldTimeSeries(filename1, "B")
Dt = FieldTimeSeries(filename1, "D")
Ft = FieldTimeSeries(filename1, "F")

N1t = FieldTimeSeries(filename3, "N")
P1t = FieldTimeSeries(filename3, "P")
D1t = FieldTimeSeries(filename3, "D")
F1t = FieldTimeSeries(filename3, "F")

t = Nt.times[1901:2500]
nt = length(t) 
z_total = znodes(Nt)
z = z_total[51:100] # only upper 500 m

N_time = zeros(nt,length(z))
P_time = zeros(nt,length(z))
B_time = zeros(nt,length(z))
D_time = zeros(nt,length(z))
F_time = zeros(nt,length(z))

N1_time = zeros(nt,length(z))
P1_time = zeros(nt,length(z))
D1_time = zeros(nt,length(z))
F1_time = zeros(nt,length(z))

for i in 1:nt
    N_time[i,:] = interior(Nt[i+1900], 1, 1, 51:100)
    P_time[i,:] = interior(Pt[i+1900], 1, 1, 51:100)
    B_time[i,:] = interior(Bt[i+1900], 1, 1, 51:100)
    D_time[i,:] = interior(Dt[i+1900], 1, 1, 51:100)
    F_time[i,:] = interior(Ft[i+1900], 1, 1, 51:100)

    N1_time[i,:] = interior(N1t[i+1900], 1, 1, 51:100)
    P1_time[i,:] = interior(P1t[i+1900], 1, 1, 51:100)
    D1_time[i,:] = interior(D1t[i+1900], 1, 1, 51:100)
    F1_time[i,:] = interior(F1t[i+1900], 1, 1, 51:100)
end

#= Start plotting
fig = Figure(;size=(1200, 600))

# First, NPBD model
# Create an Axis object
axN = Axis(fig[2, 1], xlabel = "Time (days)", ylabel = "Depth (m)",title = "[Nutrient] (mmol N m⁻³)",
            xticks = (100:200:500, string.((2000:200:2400))),yticks = (0:10:50, string.((500:-100:0))))
axP = Axis(fig[2, 3], xlabel = "Time (days)", ylabel = "Depth (m)",title = "[Phytoplankton] (mmol N m⁻³)",
            xticks = (100:200:500, string.((2000:200:2400))),yticks = (0:10:50, string.((500:-100:0))))
axD = Axis(fig[2, 5], xlabel = "Time (days)", ylabel = "Depth (m)",title = "[Detritus] (mmol N m⁻³)",
            xticks = (100:200:500, string.((2000:200:2400))),yticks = (0:10:50, string.((500:-100:0))))
axB = Axis(fig[2, 7], xlabel = "Time (days)", ylabel = "Depth (m)",title = "[Bacteria] (mmol N m⁻³)",
            xticks = (100:200:500, string.((2000:200:2400))),yticks = (0:10:50, string.((500:-100:0))))
# axF = Axis(fig[2, 9])

# Plot the heatmap
hmN = heatmap!(axN, N_time, colormap=:plasma, colorrange = (0, 40))
Colorbar(fig[2, 2], hmN)

hmP = heatmap!(axP, P_time, colormap=:plasma, colorrange = (0, 1.0))
Colorbar(fig[2, 4], hmP)

hmD = heatmap!(axD, D_time, colormap=:plasma, colorrange = (0, 2.0))
Colorbar(fig[2, 6], hmD)

hmB = heatmap!(axB, B_time, colormap=:plasma, colorrange = (0, 0.2))
Colorbar(fig[2, 8], hmB)

# hmF = heatmap!(axF, F_time, colormap=:plasma, colorrange = (0, 0.007))
# axF.xlabel = "Time"
# axF.ylabel = "Depth"
# axF.title = "[Fe] (mmol m⁻³)"
# Colorbar(fig[2, 10], hmF)

# Then, NPD model
axN1 = Axis(fig[1, 1], xlabel = "Time (days)", ylabel = "Depth (m)",title = "[Nutrient] (mmol N m⁻³)",
            xticks = (100:200:500, string.((2000:200:2400))),yticks = (0:10:50, string.((500:-100:0))))
axP1 = Axis(fig[1, 3], xlabel = "Time (days)", ylabel = "Depth (m)",title = "[Phytoplankton] (mmol N m⁻³)",
            xticks = (100:200:500, string.((2000:200:2400))),yticks = (0:10:50, string.((500:-100:0))))
axD1 = Axis(fig[1, 5], xlabel = "Time (days)", ylabel = "Depth (m)",title = "[Detritus] (mmol N m⁻³)",
            xticks = (100:200:500, string.((2000:200:2400))),yticks = (0:10:50, string.((500:-100:0))))
# axF1 = Axis(fig[1, 9])

hmN1 = heatmap!(axN1, N1_time, colormap=:plasma, colorrange = (0, 40))
Colorbar(fig[1, 2], hmN1)

hmP1 = heatmap!(axP1, P1_time, colormap=:plasma, colorrange = (0, 1.0))
Colorbar(fig[1, 4], hmP1)

hmD1 = heatmap!(axD1, D1_time, colormap=:plasma, colorrange = (0, 2.0))
Colorbar(fig[1, 6], hmD1)

# hmF1 = heatmap!(axF1, F1_time, colormap=:plasma, colorrange = (0, 0.007))
# axF1.xlabel = "Time"
# axF1.ylabel = "Depth"
# axF1.title = "[Fe] (mmol m⁻³)"
# Colorbar(fig[1, 10], hmF1)

# Display the figure
# display(fig)

# Save the last frame as a figure
# save("test.png", fig)
=#

remin_NPBD = 0.9.*B_time.*min.(D_time./(D_time.+0.05),F_time./(F_time.+1e-4)) # day-1

remin_NPD = 0.1.*D1_time

fig2 = Figure(;size=(700, 400))
axr = Axis(fig2[1, 3],xlabel = "Time (days)", ylabel = "Depth (m)",title = "NPBD remin. rate (mmol N m⁻³ day⁻¹)",
            xticks = (100:200:500, string.((2000:200:2400))),yticks = (0:10:50, string.((500:-100:0))))

hmr = heatmap!(axr, remin_NPBD, colormap=:matter, colorrange = (0, 0.25))
Colorbar(fig2[1, 4], hmr)

axr1 = Axis(fig2[1, 1],xlabel = "Time (days)", ylabel = "Depth (m)",title = "NPD remin. rate (mmol N m⁻³ day⁻¹)",
            xticks = (100:200:500, string.((2000:200:2400))),yticks = (0:10:50, string.((500:-100:0))))

hmr1 = heatmap!(axr1, remin_NPD, colormap=:matter, colorrange = (0, 0.25))
Colorbar(fig2[1, 2], hmr1)

display(fig2)
save("test2.png", fig2)