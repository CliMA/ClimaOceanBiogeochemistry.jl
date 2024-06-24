# Plotting codes: visualize 1-D column model results

# Load data
using Oceananigans
using CairoMakie

filename1 = "po4_dop_pop_sinkdown1.jld2"
filename2 = "po4_dop_pop_sinkup1.jld2"

POP1 = FieldTimeSeries(filename1, "POP")
POP2 = FieldTimeSeries(filename2, "POP")

t = POP1.times
nt = length(t) 
z = znodes(POP1)

P1 = zeros(nt,length(z))
P2 = zeros(nt,length(z))

for i in 1:nt
    P_ave[i,:] = (interior(P1[i], 1, 1, :).+interior(P2[i], 1, 1, :))./2
end

# Start plotting
fig = Figure(;size=(1200, 600))

# First, NPBD model
# Create an Axis object
axP = Axis(fig[1, 1])

# Plot the heatmap
lines!(axP, P_ave, z)
axP.xlabel = "[POP] (mmol m⁻³)"
axN.ylabel = "Depth (m)"

# Display the figure
display(fig)

# Save the last frame as a figure
save("test.png", fig)