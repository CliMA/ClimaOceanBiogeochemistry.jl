# Plotting codes: visualize 1-D column model results
# NPD vs NPBD models

# Load data
using Oceananigans
using GLMakie

filename_ctrl1 = "plot/NPBD_3.jld2"
filename_ctrl2 = "plot/NPD_2.jld2"

Nt_ctrl = FieldTimeSeries(filename_ctrl1, "N")
Pt_ctrl = FieldTimeSeries(filename_ctrl1, "P")
Bt_ctrl = FieldTimeSeries(filename_ctrl1, "B")
Dt_ctrl = FieldTimeSeries(filename_ctrl1, "D")

N1t_ctrl = FieldTimeSeries(filename_ctrl2, "N")
P1t_ctrl = FieldTimeSeries(filename_ctrl2, "P")
D1t_ctrl = FieldTimeSeries(filename_ctrl2, "D")

filename1 = "plot/pert_NPBD_130m_M05.jld2"
filename2 = "plot/pert_NPD_130m_M05.jld2"

Nt = FieldTimeSeries(filename1, "N")
Pt = FieldTimeSeries(filename1, "P")
Bt = FieldTimeSeries(filename1, "B")
Dt = FieldTimeSeries(filename1, "D")

N1t = FieldTimeSeries(filename2, "N")
P1t = FieldTimeSeries(filename2, "P")
D1t = FieldTimeSeries(filename2, "D")

filename3 = "plot/pert_NPBD_130m_P05.jld2"
filename4 = "plot/pert_NPD_130m_P05.jld2"

N2t = FieldTimeSeries(filename3, "N")
P2t = FieldTimeSeries(filename3, "P")
B2t = FieldTimeSeries(filename3, "B")
D2t = FieldTimeSeries(filename3, "D")

N3t = FieldTimeSeries(filename4, "N")
P3t = FieldTimeSeries(filename4, "P")
D3t = FieldTimeSeries(filename4, "D")

t = Nt.times
nt = length(t)
z = znodes(Nt)

N_diff = zeros(length(t),length(z))
P_diff = zeros(length(t),length(z))
B_diff = zeros(length(t),length(z))
D_diff = zeros(length(t),length(z))

N1_diff = zeros(length(t),length(z))
P1_diff = zeros(length(t),length(z))
D1_diff = zeros(length(t),length(z))

for i in 1:nt
    N_diff[i,:] = (interior(Nt[i], 1, 1, :)+interior(N2t[i], 1, 1, :))/2-interior(Nt_ctrl[i], 1, 1, :)
    P_diff[i,:] = (interior(Pt[i], 1, 1, :)+interior(P2t[i], 1, 1, :))/2-interior(Pt_ctrl[i], 1, 1, :)
    B_diff[i,:] = (interior(Bt[i], 1, 1, :)+interior(B2t[i], 1, 1, :))/2-interior(Bt_ctrl[i], 1, 1, :)
    D_diff[i,:] = (interior(Dt[i], 1, 1, :)+interior(D2t[i], 1, 1, :))/2-interior(Dt_ctrl[i], 1, 1, :)

    N1_diff[i,:] = (interior(N1t[i], 1, 1, :)+interior(N3t[i], 1, 1, :))/2-interior(N1t_ctrl[i], 1, 1, :)
    P1_diff[i,:] = (interior(P1t[i], 1, 1, :)+interior(P3t[i], 1, 1, :))/2-interior(P1t_ctrl[i], 1, 1, :)
    D1_diff[i,:] = (interior(D1t[i], 1, 1, :)+interior(D3t[i], 1, 1, :))/2-interior(D1t_ctrl[i], 1, 1, :)
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
hmN = heatmap!(axN, N_diff, colormap=:bwr, colorrange = (-0.1, 0.1))
axN.xlabel = "Time"
axN.ylabel = "Depth"
axN.title = "[Nutrient] (mmol m⁻³)"
Colorbar(fig[1, 2], hmN)

hmP = heatmap!(axP, P_diff, colormap=:bwr, colorrange = (-0.01, 0.01))
axP.xlabel = "Time"
axP.ylabel = "Depth"
axP.title = "[Phytoplankton] (mmol m⁻³)"
Colorbar(fig[1, 4], hmP)

hmB = heatmap!(axB, B_diff, colormap=:bwr, colorrange = (-0.01, 0.01))
axB.xlabel = "Time"
axB.ylabel = "Depth"
axB.title = "[Bacteria] (mmol m⁻³)"
Colorbar(fig[1, 6], hmB)

hmD = heatmap!(axD, D_diff, colormap=:bwr, colorrange = (-0.01, 0.01))
axD.xlabel = "Time"
axD.ylabel = "Depth"
axD.title = "[Detritus] (mmol m⁻³)"
Colorbar(fig[1, 8], hmD)

# Then, NPD model
axN1 = Axis(fig[2, 1])
axP1 = Axis(fig[2, 3])
axD1 = Axis(fig[2, 7])

hmN1 = heatmap!(axN1, N1_diff, colormap=:bwr, colorrange = (-0.1, 0.1))
axN1.xlabel = "Time"
axN1.ylabel = "Depth"
axN1.title = "[Nutrient] (mmol m⁻³)"
Colorbar(fig[2, 2], hmN1)

hmP1 = heatmap!(axP1, P1_diff, colormap=:bwr, colorrange = (-0.01, 0.01))
axP1.xlabel = "Time"
axP1.ylabel = "Depth"
axP1.title = "[Phytoplankton] (mmol m⁻³)"
Colorbar(fig[2, 4], hmP1)

hmD1 = heatmap!(axD1, D1_diff, colormap=:bwr, colorrange = (-0.01, 0.01))
axD1.xlabel = "Time"
axD1.ylabel = "Depth"
axD1.title = "[Detritus] (mmol m⁻³)"
Colorbar(fig[2, 8], hmD1)

# Display the figure
display(fig)

# Save the last frame as a figure
save("test.png", fig)