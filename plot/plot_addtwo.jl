# Plotting codes: visualize 1-D column model results
# NPD vs NPBD models

# Load data
using Oceananigans
using GLMakie

filename_ctrl1 = "NPBD_4.jld2"
filename_ctrl2 = "NPD_3.jld2"

Nt_ctrl = FieldTimeSeries(filename_ctrl1, "N")
Pt_ctrl = FieldTimeSeries(filename_ctrl1, "P")
Bt_ctrl = FieldTimeSeries(filename_ctrl1, "B")
Dt_ctrl = FieldTimeSeries(filename_ctrl1, "D")

N1t_ctrl = FieldTimeSeries(filename_ctrl2, "N")
P1t_ctrl = FieldTimeSeries(filename_ctrl2, "P")
D1t_ctrl = FieldTimeSeries(filename_ctrl2, "D")

# shoaling
filename01 = "plot/NPBD_MLD1_before.jld2"
filename02 = "plot/NPBD_MLD1_after.jld2"

filename03 = "plot/NPD_MLD1_before.jld2"
filename04 = "plot/NPD_MLD1_after.jld2"

N0t = FieldTimeSeries(filename01, "N")
P0t = FieldTimeSeries(filename01, "P")
B0t = FieldTimeSeries(filename01, "B")
D0t = FieldTimeSeries(filename01, "D")

N0t_after = FieldTimeSeries(filename02, "N")
P0t_after = FieldTimeSeries(filename02, "P")
B0t_after = FieldTimeSeries(filename02, "B")
D0t_after = FieldTimeSeries(filename02, "D")

N01t = FieldTimeSeries(filename03, "N")
P01t = FieldTimeSeries(filename03, "P")
D01t = FieldTimeSeries(filename03, "D")

N01t_after = FieldTimeSeries(filename04, "N")
P01t_after = FieldTimeSeries(filename04, "P")
D01t_after = FieldTimeSeries(filename04, "D")

# deepening
filename1 = "plot/NPBD_MLD2_before.jld2"
filename2 = "plot/NPBD_MLD2_after.jld2"

filename3 = "plot/NPD_MLD2_before.jld2"
filename4 = "plot/NPD_MLD2_after.jld2"

Nt = FieldTimeSeries(filename1, "N")
Pt = FieldTimeSeries(filename1, "P")
Bt = FieldTimeSeries(filename1, "B")
Dt = FieldTimeSeries(filename1, "D")

Nt_after = FieldTimeSeries(filename2, "N")
Pt_after = FieldTimeSeries(filename2, "P")
Bt_after = FieldTimeSeries(filename2, "B")
Dt_after = FieldTimeSeries(filename2, "D")

N1t = FieldTimeSeries(filename3, "N")
P1t = FieldTimeSeries(filename3, "P")
D1t = FieldTimeSeries(filename3, "D")

N1t_after = FieldTimeSeries(filename4, "N")
P1t_after = FieldTimeSeries(filename4, "P")
D1t_after = FieldTimeSeries(filename4, "D")

t = Nt.times
nt = length(t) *2
z = znodes(Nt)

N0_time = zeros(nt,length(z))
P0_time = zeros(nt,length(z))
B0_time = zeros(nt,length(z))
D0_time = zeros(nt,length(z))

N01_time = zeros(nt,length(z))
P01_time = zeros(nt,length(z))
D01_time = zeros(nt,length(z))

N_time = zeros(nt,length(z))
P_time = zeros(nt,length(z))
B_time = zeros(nt,length(z))
D_time = zeros(nt,length(z))

N1_time = zeros(nt,length(z))
P1_time = zeros(nt,length(z))
D1_time = zeros(nt,length(z))

for i in 1:length(t)
    # NPBD_shoaling
    N0_time[i,:] = interior(N0t[i], 1, 1, :)
    P0_time[i,:] = interior(P0t[i], 1, 1, :)
    B0_time[i,:] = interior(B0t[i], 1, 1, :)
    D0_time[i,:] = interior(D0t[i], 1, 1, :)

    N0_time[i+length(t),:] = interior(N0t_after[i], 1, 1, :)
    P0_time[i+length(t),:] = interior(P0t_after[i], 1, 1, :)
    B0_time[i+length(t),:] = interior(B0t_after[i], 1, 1, :)
    D0_time[i+length(t),:] = interior(D0t_after[i], 1, 1, :)

    N01_time[i,:] = interior(N01t[i], 1, 1, :)
    P01_time[i,:] = interior(P01t[i], 1, 1, :)
    D01_time[i,:] = interior(D01t[i], 1, 1, :)

    N01_time[i+length(t),:] = interior(N01t_after[i], 1, 1, :)
    P01_time[i+length(t),:] = interior(P01t_after[i], 1, 1, :)
    D01_time[i+length(t),:] = interior(D01t_after[i], 1, 1, :)

        # NPBD_deepening
        N_time[i,:] = interior(Nt[i], 1, 1, :)
        P_time[i,:] = interior(Pt[i], 1, 1, :)
        B_time[i,:] = interior(Bt[i], 1, 1, :)
        D_time[i,:] = interior(Dt[i], 1, 1, :)
    
        N_time[i+length(t),:] = interior(Nt_after[i], 1, 1, :)
        P_time[i+length(t),:] = interior(Pt_after[i], 1, 1, :)
        B_time[i+length(t),:] = interior(Bt_after[i], 1, 1, :)
        D_time[i+length(t),:] = interior(Dt_after[i], 1, 1, :)
    
        N1_time[i,:] = interior(N1t[i], 1, 1, :)
        P1_time[i,:] = interior(P1t[i], 1, 1, :)
        D1_time[i,:] = interior(D1t[i], 1, 1, :)
    
        N1_time[i+length(t),:] = interior(N1t_after[i], 1, 1, :)
        P1_time[i+length(t),:] = interior(P1t_after[i], 1, 1, :)
        D1_time[i+length(t),:] = interior(D1t_after[i], 1, 1, :)
end

N_diff = zeros(nt,length(z))
P_diff = zeros(nt,length(z))
B_diff = zeros(nt,length(z))
D_diff = zeros(nt,length(z))

N1_diff = zeros(nt,length(z))
P1_diff = zeros(nt,length(z))
D1_diff = zeros(nt,length(z))

for i in 1:(nt-1)
    N_diff[i,:] = (N0_time[i,:]+N_time[i,:])/2-interior(Nt_ctrl[i], 1, 1, :)
    P_diff[i,:] = (P0_time[i,:]+P_time[i,:])/2-interior(Pt_ctrl[i], 1, 1, :)
    B_diff[i,:] = (B0_time[i,:]+B_time[i,:])/2-interior(Bt_ctrl[i], 1, 1, :)
    D_diff[i,:] = (D0_time[i,:]+D_time[i,:])/2-interior(Dt_ctrl[i], 1, 1, :)

    N1_diff[i,:] = (N01_time[i,:]+N1_time[i,:])/2-interior(N1t_ctrl[i], 1, 1, :)
    P1_diff[i,:] = (P01_time[i,:]+P1_time[i,:])/2-interior(P1t_ctrl[i], 1, 1, :)
    D1_diff[i,:] = (D01_time[i,:]+D1_time[i,:])/2-interior(D1t_ctrl[i], 1, 1, :)
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

hmP = heatmap!(axP, P_diff, colormap=:bwr, colorrange = (-0.05, 0.05))
axP.xlabel = "Time"
axP.ylabel = "Depth"
axP.title = "[Phytoplankton] (mmol m⁻³)"
Colorbar(fig[1, 4], hmP)

hmB = heatmap!(axB, B_diff, colormap=:bwr, colorrange = (-0.05, 0.05))
axB.xlabel = "Time"
axB.ylabel = "Depth"
axB.title = "[Bacteria] (mmol m⁻³)"
Colorbar(fig[1, 6], hmB)

hmD = heatmap!(axD, D_diff, colormap=:bwr, colorrange = (-0.2, 0.2))
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

hmP1 = heatmap!(axP1, P1_diff, colormap=:bwr, colorrange = (-0.05, 0.05))
axP1.xlabel = "Time"
axP1.ylabel = "Depth"
axP1.title = "[Phytoplankton] (mmol m⁻³)"
Colorbar(fig[2, 4], hmP1)

hmD1 = heatmap!(axD1, D1_diff, colormap=:bwr, colorrange = (-0.2, 0.2))
axD1.xlabel = "Time"
axD1.ylabel = "Depth"
axD1.title = "[Detritus] (mmol m⁻³)"
Colorbar(fig[2, 8], hmD1)

# Display the figure
display(fig)

# Save the last frame as a figure
save("test.png", fig)