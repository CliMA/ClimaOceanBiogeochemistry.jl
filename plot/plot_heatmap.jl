# Plotting codes: visualize 1-D column model results
# NPD vs NPBD models

# Load data
using Oceananigans
using GLMakie

filename1 = "plot/pert_NPBD_130m_P05.jld2"
filename2 = "plot/pert_NPD_130m_P05.jld2"
# filename1 = "plot/NPBD_3.jld2"
# filename2 = "plot/NPD_2.jld2"

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

N_time = zeros(length(t),length(z))
P_time = zeros(length(t),length(z))
B_time = zeros(length(t),length(z))
D_time = zeros(length(t),length(z))

N1_time = zeros(length(t),length(z))
P1_time = zeros(length(t),length(z))
D1_time = zeros(length(t),length(z))

for i in 1:nt
    N_time[i,:] = interior(Nt[i], 1, 1, :)
    P_time[i,:] = interior(Pt[i], 1, 1, :)
    B_time[i,:] = interior(Bt[i], 1, 1, :)
    D_time[i,:] = interior(Dt[i], 1, 1, :)

    N1_time[i,:] = interior(N1t[i], 1, 1, :)
    P1_time[i,:] = interior(P1t[i], 1, 1, :)
    D1_time[i,:] = interior(D1t[i], 1, 1, :)
end

# Start plotting
fig = Figure(;size=(1200, 600))

# First, NPBD model
# Create an Axis object
axN = Axis(fig[1, 1])
axP = Axis(fig[1, 3])
axB = Axis(fig[1, 5])
axD = Axis(fig[1, 7])

# Plot the heatmap
hmN = heatmap!(axN, N_time, colormap=:plasma, colorrange = (0, 30))
axN.xlabel = "Time"
axN.ylabel = "Depth"
axN.title = "[Nutrient] (mmol m⁻³)"
Colorbar(fig[1, 2], hmN)

hmP = heatmap!(axP, P_time, colormap=:plasma, colorrange = (0, 0.5))
axP.xlabel = "Time"
axP.ylabel = "Depth"
axP.title = "[Phytoplankton] (mmol m⁻³)"
Colorbar(fig[1, 4], hmP)

hmB = heatmap!(axB, B_time, colormap=:plasma, colorrange = (0, 0.2))
axB.xlabel = "Time"
axB.ylabel = "Depth"
axB.title = "[Bacteria] (mmol m⁻³)"
Colorbar(fig[1, 6], hmB)

hmD = heatmap!(axD, D_time, colormap=:plasma, colorrange = (0, 1))
axD.xlabel = "Time"
axD.ylabel = "Depth"
axD.title = "[Detritus] (mmol m⁻³)"
Colorbar(fig[1, 8], hmD)

# Then, NPD model
axN1 = Axis(fig[2, 1])
axP1 = Axis(fig[2, 3])
axD1 = Axis(fig[2, 7])

hmN1 = heatmap!(axN1, N1_time, colormap=:plasma, colorrange = (0, 30))
axN1.xlabel = "Time"
axN1.ylabel = "Depth"
axN1.title = "[Nutrient] (mmol m⁻³)"
Colorbar(fig[2, 2], hmN1)

hmP1 = heatmap!(axP1, P1_time, colormap=:plasma, colorrange = (0, 0.5))
axP1.xlabel = "Time"
axP1.ylabel = "Depth"
axP1.title = "[Phytoplankton] (mmol m⁻³)"
Colorbar(fig[2, 4], hmP1)

hmD1 = heatmap!(axD1, D1_time, colormap=:plasma, colorrange = (0, 1))
axD1.xlabel = "Time"
axD1.ylabel = "Depth"
axD1.title = "[Detritus] (mmol m⁻³)"
Colorbar(fig[2, 8], hmD1)

# Display the figure
display(fig)

# Save the last frame as a figure
save("test.png", fig)