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

t = Nt.times[1001:3000]
nt = length(t) 
z = znodes(Nt)

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
    N_time[i,:] = interior(Nt[i+1000], 1, 1, :)
    P_time[i,:] = interior(Pt[i+1000], 1, 1, :)
    B_time[i,:] = interior(Bt[i+1000], 1, 1, :)
    D_time[i,:] = interior(Dt[i+1000], 1, 1, :)
    F_time[i,:] = interior(Ft[i+1000], 1, 1, :)

    N1_time[i,:] = interior(N1t[i+1000], 1, 1, :)
    P1_time[i,:] = interior(P1t[i+1000], 1, 1, :)
    D1_time[i,:] = interior(D1t[i+1000], 1, 1, :)
    F1_time[i,:] = interior(F1t[i+1000], 1, 1, :)
end

# Start plotting
fig = Figure(;size=(1300, 600))

# First, NPBD model
# Create an Axis object
axN = Axis(fig[2, 1])
axP = Axis(fig[2, 3])
axD = Axis(fig[2, 5])
axB = Axis(fig[2, 7])
axF = Axis(fig[2, 9])

# Plot the heatmap
hmN = heatmap!(axN, N_time, colormap=:plasma, colorrange = (0, 40))
axN.xlabel = "Time"
axN.ylabel = "Depth"
axN.title = "[Nutrient] (mmol m⁻³)"
Colorbar(fig[2, 2], hmN)

hmP = heatmap!(axP, P_time, colormap=:plasma, colorrange = (0, 1.0))
axP.xlabel = "Time"
axP.ylabel = "Depth"
axP.title = "[Phytoplankton] (mmol m⁻³)"
Colorbar(fig[2, 4], hmP)

hmD = heatmap!(axD, D_time, colormap=:plasma, colorrange = (0, 2.0))
axD.xlabel = "Time"
axD.ylabel = "Depth"
axD.title = "[Detritus] (mmol m⁻³)"
Colorbar(fig[2, 6], hmD)

hmB = heatmap!(axB, B_time, colormap=:plasma, colorrange = (0, 0.2))
axB.xlabel = "Time"
axB.ylabel = "Depth"
axB.title = "[Bacteria] (mmol m⁻³)"
Colorbar(fig[2, 8], hmB)

hmF = heatmap!(axF, F_time, colormap=:plasma, colorrange = (0, 0.007))
axF.xlabel = "Time"
axF.ylabel = "Depth"
axF.title = "[Fe] (mmol m⁻³)"
Colorbar(fig[2, 10], hmF)

# Then, NPD model
axN1 = Axis(fig[1, 1])
axP1 = Axis(fig[1, 3])
axD1 = Axis(fig[1, 5])
axF1 = Axis(fig[1, 9])

hmN1 = heatmap!(axN1, N1_time, colormap=:plasma, colorrange = (0, 40))
axN1.xlabel = "Time"
axN1.ylabel = "Depth"
axN1.title = "[Nutrient] (mmol m⁻³)"
Colorbar(fig[1, 2], hmN1)

hmP1 = heatmap!(axP1, P1_time, colormap=:plasma, colorrange = (0, 1.0))
axP1.xlabel = "Time"
axP1.ylabel = "Depth"
axP1.title = "[Phytoplankton] (mmol m⁻³)"
Colorbar(fig[1, 4], hmP1)

hmD1 = heatmap!(axD1, D1_time, colormap=:plasma, colorrange = (0, 2.0))
axD1.xlabel = "Time"
axD1.ylabel = "Depth"
axD1.title = "[Detritus] (mmol m⁻³)"
Colorbar(fig[1, 6], hmD1)

hmF1 = heatmap!(axF1, F1_time, colormap=:plasma, colorrange = (0, 0.007))
axF1.xlabel = "Time"
axF1.ylabel = "Depth"
axF1.title = "[Fe] (mmol m⁻³)"
Colorbar(fig[1, 10], hmF1)

# Display the figure
display(fig)

# Save the last frame as a figure
save("test.png", fig)