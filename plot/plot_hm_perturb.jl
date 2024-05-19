# Plotting codes: visualize 1-D column model results
# NPD vs NPBD models

# Load data
using Oceananigans
using GLMakie

filename1 = "plot/NPBD_FeFert1.jld2"
filename2 = "plot/NPBD_FeFert2.jld2"
filename3 = "plot/NPBD_FeFert3.jld2"

filename4 = "plot/NPD_FeFert1.jld2"
filename5 = "plot/NPD_FeFert2.jld2"
filename6 = "plot/NPD_FeFert3.jld2"

Nt = FieldTimeSeries(filename1, "N")
Pt = FieldTimeSeries(filename1, "P")
Bt = FieldTimeSeries(filename1, "B")
Dt = FieldTimeSeries(filename1, "D")

Nt_2 = FieldTimeSeries(filename2, "N")
Pt_2 = FieldTimeSeries(filename2, "P")
Bt_2 = FieldTimeSeries(filename2, "B")
Dt_2 = FieldTimeSeries(filename2, "D")

Nt_3 = FieldTimeSeries(filename3, "N")
Pt_3 = FieldTimeSeries(filename3, "P")
Bt_3 = FieldTimeSeries(filename3, "B")
Dt_3 = FieldTimeSeries(filename3, "D")

Nt_4 = FieldTimeSeries(filename4, "N")
Pt_4 = FieldTimeSeries(filename4, "P")
Dt_4 = FieldTimeSeries(filename4, "D")

Nt_5 = FieldTimeSeries(filename5, "N")
Pt_5 = FieldTimeSeries(filename5, "P")
Dt_5 = FieldTimeSeries(filename5, "D")

Nt_6 = FieldTimeSeries(filename6, "N")
Pt_6 = FieldTimeSeries(filename6, "P")
Dt_6 = FieldTimeSeries(filename6, "D")

# nt = length(Nt.times)+ length(Nt_2.times)+ length(Nt_3.times)
nt = 300 + length(Nt_2.times)+300 
z = znodes(Nt)

N_time = zeros(nt,length(z))
P_time = zeros(nt,length(z))
B_time = zeros(nt,length(z))
D_time = zeros(nt,length(z))

N1_time = zeros(nt,length(z))
P1_time = zeros(nt,length(z))
D1_time = zeros(nt,length(z))

# for i in 1:length(Nt.times)
#     N_time[i,:] = interior(Nt[i], 1, 1, :)
#     P_time[i,:] = interior(Pt[i], 1, 1, :)
#     B_time[i,:] = interior(Bt[i], 1, 1, :)
#     D_time[i,:] = interior(Dt[i], 1, 1, :)

#     N1_time[i,:] = interior(Nt_4[i], 1, 1, :)
#     P1_time[i,:] = interior(Pt_4[i], 1, 1, :)
#     D1_time[i,:] = interior(Dt_4[i], 1, 1, :)
# end
# for i in 1:length(Nt_2.times)
#     N_time[length(Nt.times)+i,:] = interior(Nt_2[i], 1, 1, :)
#     P_time[length(Nt.times)+i,:] = interior(Pt_2[i], 1, 1, :)
#     B_time[length(Nt.times)+i,:] = interior(Bt_2[i], 1, 1, :)
#     D_time[length(Nt.times)+i,:] = interior(Dt_2[i], 1, 1, :)

#     N1_time[length(Nt.times)+i,:] = interior(Nt_5[i], 1, 1, :)
#     P1_time[length(Nt.times)+i,:] = interior(Pt_5[i], 1, 1, :)
#     D1_time[length(Nt.times)+i,:] = interior(Dt_5[i], 1, 1, :)
# end
# for i in 1:length(Nt_3.times)
#     N_time[length(Nt.times)+length(Nt_2.times)+i,:] = interior(Nt_3[i], 1, 1, :)
#     P_time[length(Nt.times)+length(Nt_2.times)+i,:] = interior(Pt_3[i], 1, 1, :)
#     B_time[length(Nt.times)+length(Nt_2.times)+i,:] = interior(Bt_3[i], 1, 1, :)
#     D_time[length(Nt.times)+length(Nt_2.times)+i,:] = interior(Dt_3[i], 1, 1, :)

#     N1_time[length(Nt.times)+length(Nt_2.times)+i,:] = interior(Nt_6[i], 1, 1, :)
#     P1_time[length(Nt.times)+length(Nt_2.times)+i,:] = interior(Pt_6[i], 1, 1, :)
#     D1_time[length(Nt.times)+length(Nt_2.times)+i,:] = interior(Dt_6[i], 1, 1, :)
# end
for i in 1:300
    N_time[i,:] = interior(Nt[i], 1, 1, :)
    P_time[i,:] = interior(Pt[i], 1, 1, :)
    B_time[i,:] = interior(Bt[i], 1, 1, :)
    D_time[i,:] = interior(Dt[i], 1, 1, :)

    N1_time[i,:] = interior(Nt_4[i], 1, 1, :)
    P1_time[i,:] = interior(Pt_4[i], 1, 1, :)
    D1_time[i,:] = interior(Dt_4[i], 1, 1, :)
end
for i in 1:length(Nt_2.times)
    N_time[300+i,:] = interior(Nt_2[i], 1, 1, :)
    P_time[300+i,:] = interior(Pt_2[i], 1, 1, :)
    B_time[300+i,:] = interior(Bt_2[i], 1, 1, :)
    D_time[300+i,:] = interior(Dt_2[i], 1, 1, :)

    N1_time[300+i,:] = interior(Nt_5[i], 1, 1, :)
    P1_time[300+i,:] = interior(Pt_5[i], 1, 1, :)
    D1_time[300+i,:] = interior(Dt_5[i], 1, 1, :)
end
for i in 1:300
    N_time[300+length(Nt_2.times)+i,:] = interior(Nt_3[i], 1, 1, :)
    P_time[300+length(Nt_2.times)+i,:] = interior(Pt_3[i], 1, 1, :)
    B_time[300+length(Nt_2.times)+i,:] = interior(Bt_3[i], 1, 1, :)
    D_time[300+length(Nt_2.times)+i,:] = interior(Dt_3[i], 1, 1, :)

    N1_time[300+length(Nt_2.times)+i,:] = interior(Nt_6[i], 1, 1, :)
    P1_time[300+length(Nt_2.times)+i,:] = interior(Pt_6[i], 1, 1, :)
    D1_time[300+length(Nt_2.times)+i,:] = interior(Dt_6[i], 1, 1, :)
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
hmN = heatmap!(axN, N_time, colormap=:plasma, colorrange = (0, 35))
axN.xlabel = "Time"
axN.ylabel = "Depth"
axN.title = "[Nutrient] (mmol m⁻³)"
Colorbar(fig[1, 2], hmN)

hmP = heatmap!(axP, P_time, colormap=:plasma, colorrange = (0, 0.75))
axP.xlabel = "Time"
axP.ylabel = "Depth"
axP.title = "[Phytoplankton] (mmol m⁻³)"
Colorbar(fig[1, 4], hmP)

hmB = heatmap!(axB, B_time, colormap=:plasma, colorrange = (0, 0.15))
axB.xlabel = "Time"
axB.ylabel = "Depth"
axB.title = "[Bacteria] (mmol m⁻³)"
Colorbar(fig[1, 6], hmB)

hmD = heatmap!(axD, D_time, colormap=:plasma, colorrange = (0, 1.5))
axD.xlabel = "Time"
axD.ylabel = "Depth"
axD.title = "[Detritus] (mmol m⁻³)"
Colorbar(fig[1, 8], hmD)

# Then, NPD model
axN1 = Axis(fig[2, 1])
axP1 = Axis(fig[2, 3])
axD1 = Axis(fig[2, 7])

hmN1 = heatmap!(axN1, N1_time, colormap=:plasma, colorrange = (0, 35))
axN1.xlabel = "Time"
axN1.ylabel = "Depth"
axN1.title = "[Nutrient] (mmol m⁻³)"
Colorbar(fig[2, 2], hmN1)

hmP1 = heatmap!(axP1, P1_time, colormap=:plasma, colorrange = (0, 0.75))
axP1.xlabel = "Time"
axP1.ylabel = "Depth"
axP1.title = "[Phytoplankton] (mmol m⁻³)"
Colorbar(fig[2, 4], hmP1)

hmD1 = heatmap!(axD1, D1_time, colormap=:plasma, colorrange = (0, 1.5))
axD1.xlabel = "Time"
axD1.ylabel = "Depth"
axD1.title = "[Detritus] (mmol m⁻³)"
Colorbar(fig[2, 8], hmD1)

# Display the figure
display(fig)

# Save the last frame as a figure
save("test.png", fig)