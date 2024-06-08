# Plotting codes: visualize 1-D column model results
# NPD vs NPBD models

# Load data
using Oceananigans
using GLMakie

filename_ctrl1 = "plot/NPBD_5.jld2"
filename_ctrl2 = "plot/NPD_4.jld2"

Nt_ctrl = FieldTimeSeries(filename_ctrl1, "N")
Pt_ctrl = FieldTimeSeries(filename_ctrl1, "P")
Bt_ctrl = FieldTimeSeries(filename_ctrl1, "B")
Dt_ctrl = FieldTimeSeries(filename_ctrl1, "D")

N1t_ctrl = FieldTimeSeries(filename_ctrl2, "N")
P1t_ctrl = FieldTimeSeries(filename_ctrl2, "P")
D1t_ctrl = FieldTimeSeries(filename_ctrl2, "D")

# shoaling
filename01 = "plot/NPBD5_V15.jld2"
filename02 = "plot/NPBD5_V09.jld2"

filename03 = "plot/NPD4_r018.jld2"
filename04 = "plot/NPD4_r012.jld2"

# Increase
N0t = FieldTimeSeries(filename01, "N")
P0t = FieldTimeSeries(filename01, "P")
B0t = FieldTimeSeries(filename01, "B")
D0t = FieldTimeSeries(filename01, "D")
# Decrease
N0t_after = FieldTimeSeries(filename02, "N")
P0t_after = FieldTimeSeries(filename02, "P")
B0t_after = FieldTimeSeries(filename02, "B")
D0t_after = FieldTimeSeries(filename02, "D")

# Increase
N01t = FieldTimeSeries(filename03, "N")
P01t = FieldTimeSeries(filename03, "P")
D01t = FieldTimeSeries(filename03, "D")
# Decrease
N01t_after = FieldTimeSeries(filename04, "N")
P01t_after = FieldTimeSeries(filename04, "P")
D01t_after = FieldTimeSeries(filename04, "D")

t = Nt_ctrl.times
nt = length(t) 
z = znodes(Nt_ctrl)

N_diff = zeros(nt,length(z))
P_diff = zeros(nt,length(z))
B_diff = zeros(nt,length(z))
D_diff = zeros(nt,length(z))

N1_diff = zeros(nt,length(z))
P1_diff = zeros(nt,length(z))
D1_diff = zeros(nt,length(z))

for i in 1:nt
    N_diff[i,:] = (interior(N0t[i], 1, 1, :)+interior(N0t_after[i], 1, 1, :))/2-interior(Nt_ctrl[i], 1, 1, :)
    P_diff[i,:] = (interior(P0t[i], 1, 1, :)+interior(P0t_after[i], 1, 1, :))/2-interior(Pt_ctrl[i], 1, 1, :)
    B_diff[i,:] = (interior(B0t[i], 1, 1, :)+interior(B0t_after[i], 1, 1, :))/2-interior(Bt_ctrl[i], 1, 1, :)
    D_diff[i,:] = (interior(D0t[i], 1, 1, :)+interior(D0t_after[i], 1, 1, :))/2-interior(Dt_ctrl[i], 1, 1, :)

    N1_diff[i,:] = (interior(N01t[i], 1, 1, :)+interior(N01t_after[i], 1, 1, :))/2-interior(N1t_ctrl[i], 1, 1, :)
    P1_diff[i,:] = (interior(P01t[i], 1, 1, :)+interior(P01t_after[i], 1, 1, :))/2-interior(P1t_ctrl[i], 1, 1, :)
    D1_diff[i,:] = (interior(D01t[i], 1, 1, :)+interior(D01t_after[i], 1, 1, :))/2-interior(D1t_ctrl[i], 1, 1, :)
end

# Start plotting
fig = Figure(;size=(1200, 600))

# First, NPBD model
# Create an Axis object
axN = Axis(fig[2, 1])
axP = Axis(fig[2, 3])
axD = Axis(fig[2, 5])
axB = Axis(fig[2, 7])

# Plot the heatmap
hmN = heatmap!(axN, N_diff, colormap=:bwr, colorrange = (-0.1, 0.1))
axN.xlabel = "Time"
axN.ylabel = "Depth"
axN.title = "[Nutrient] (mmol m⁻³)"
Colorbar(fig[2, 2], hmN)

hmP = heatmap!(axP, P_diff, colormap=:bwr, colorrange = (-0.05, 0.05))
axP.xlabel = "Time"
axP.ylabel = "Depth"
axP.title = "[Phytoplankton] (mmol m⁻³)"
Colorbar(fig[2, 4], hmP)

hmB = heatmap!(axB, B_diff, colormap=:bwr, colorrange = (-0.05, 0.05))
axB.xlabel = "Time"
axB.ylabel = "Depth"
axB.title = "[Bacteria] (mmol m⁻³)"
Colorbar(fig[2, 8], hmB)

hmD = heatmap!(axD, D_diff, colormap=:bwr, colorrange = (-0.2, 0.2))
axD.xlabel = "Time"
axD.ylabel = "Depth"
axD.title = "[Detritus] (mmol m⁻³)"
Colorbar(fig[2, 6], hmD)

# Then, NPD model
axN1 = Axis(fig[1, 1])
axP1 = Axis(fig[1, 3])
axD1 = Axis(fig[1, 5])

hmN1 = heatmap!(axN1, N1_diff, colormap=:bwr, colorrange = (-0.1, 0.1))
axN1.xlabel = "Time"
axN1.ylabel = "Depth"
axN1.title = "[Nutrient] (mmol m⁻³)"
Colorbar(fig[1, 2], hmN1)

hmP1 = heatmap!(axP1, P1_diff, colormap=:bwr, colorrange = (-0.05, 0.05))
axP1.xlabel = "Time"
axP1.ylabel = "Depth"
axP1.title = "[Phytoplankton] (mmol m⁻³)"
Colorbar(fig[1, 4], hmP1)

hmD1 = heatmap!(axD1, D1_diff, colormap=:bwr, colorrange = (-0.2, 0.2))
axD1.xlabel = "Time"
axD1.ylabel = "Depth"
axD1.title = "[Detritus] (mmol m⁻³)"
Colorbar(fig[1, 6], hmD1)

# Display the figure
display(fig)

# Save the last frame as a figure
save("test.png", fig)