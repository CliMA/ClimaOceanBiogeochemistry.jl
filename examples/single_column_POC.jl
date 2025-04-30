# Build a one-dimensional single column model, 
# simulating POC remineralization and sinking in a minimalistic manner
# Purpose: Test the remineralization formulation of POC as a function of depth

using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Fields: ConstantField, ZeroField
using Oceananigans.Grids: Center, znode
using GLMakie
using Printf

######################################### Grid ##########################################
# the 1D domain
H = 1000 # deepest depth 
Nz = 200 # number of vertical grids (resolution = H/Nz)

grid = RectilinearGrid(size = Nz; z = (-H, 0), topology = (Flat, Flat, Bounded))

################################### Boundary condition ###################################
# Add a top boundary condition (fixed value) to ensure there is POC production from the top 

#top_value = 1.0
top_value(t) = (1.0 +sinpi(2 / 365 * (t/day))) 
value_top_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(top_value))

####################################### Forcing #######################################
# Particle forcing: remineralization and sinking

# Calculate remineralization of particulate organic phosphorus according to a first-order rate constant.
z₀ = -110 # Reference depth (m)
wₛ = -10/day # Sinking velocity (s⁻¹)
martin_b = 0.84

sinking = AdvectiveForcing(w=wₛ)

remin_func(i, j, k, grid, clock, model_fields, parameters) = -parameters.b * parameters.wₛ / (grid.zᵃᵃᶜ[k] + parameters.z₀) * model_fields.POC[i,j,k]

remineralization = Forcing(remin_func, parameters=(b = martin_b, wₛ = wₛ, z₀ = z₀), discrete_form=true) 

####################################### Model #######################################

model = HydrostaticFreeSurfaceModel(; grid,
                                    boundary_conditions = (; POC = value_top_bcs),
                                    velocities = nothing,
                                    tracers = (:POC), 
                                    forcing = (; POC = (sinking,remineralization)), 
                                    buoyancy = nothing) 

############################## Initial conditions ################################
#Dᵢ(x,y,z) = 5 * (-z).^-0.84
# POCᵢ(z) = ((z + z₀) / z₀)^(-0.84) 
set!(model, POC=1)
simulation = Simulation(model, Δt=1hour, stop_time=365*2days)

############################### Output file ###############################
filename = "POC_1D_yearly.jld2"
simulation.output_writers[:fields] = JLD2OutputWriter(model, model.tracers;
                                                      filename,
                                                      schedule = TimeInterval(1day),
                                                      overwrite_existing = true)

# Start simulation!
run!(simulation)

top_value2(t) = (1.0 +sinpi(2 / (365/12) * (t/day))) 
value_top_bcs2 = FieldBoundaryConditions(top = ValueBoundaryCondition(top_value2))

model2 = HydrostaticFreeSurfaceModel(; grid,
                                    boundary_conditions = (; POC = value_top_bcs2),
                                    velocities = nothing,
                                    tracers = (:POC), 
                                    forcing = (; POC = (sinking,remineralization)), 
                                    buoyancy = nothing) 
set!(model2, POC=1)
simulation2 = Simulation(model2, Δt=1hour, stop_time=365*2days)
filename2 = "POC_1D_monthly.jld2"
simulation2.output_writers[:fields] = JLD2OutputWriter(model2, model2.tracers;
                                                      filename = filename2,
                                                      schedule = TimeInterval(1day),
                                                      overwrite_existing = true)
run!(simulation2)

# ###################### Visualization ######################

# All that's left is to visualize the results.

POCt = FieldTimeSeries(filename, "POC")
POCt2 = FieldTimeSeries(filename2, "POC")

times = POCt.times
nt = length(times)
zw = znodes(POCt)

# Analytical solution
function D(z, t; N=365, w=-10, z0=z₀)
    return (1+sinpi(2 / N * (t - z/w))) * ((z + z0) / z0)^(-0.84) 
end
Ana_1 = [D(z, t; N=365) for z in zw, t in times/day]
Ana_2 = [D(z, t; N=365/12) for z in zw, t in times/day]

fig = Figure(;size=(800, 600))

slider = Slider(fig[2, 1:2], range=1:nt, startvalue=1)
n = slider.value
title = @lift @sprintf("t = %d days", (times[$n]/ day))
Label(fig[0, 1:2], title)
POCn = @lift interior(POCt[$n], 1, 1, :)
POCn2 = @lift interior(POCt2[$n], 1, 1, :)
POC_data = POCt[1,1,:,:]
POC_data2 = POCt2[1,1,:,:]
diff_1 = POC_data .- Ana_1 
diff_2 = POC_data2 .- Ana_2 

############################################################
######################## Yearly ###########################
############################################################

ax1 = Axis(fig[1, 1], ylabel="z (m)", xlabel="POC conc (mol m⁻³)")
xlims!(ax1, -0.5, 2.5)
ylims!(ax1, -1000, 0)
lines!(ax1, POCn, zw, linewidth = 3,label = "Modeled")
POC_prof = lines!(ax1, Ana_1[:, 1], zw, linewidth = 1.5,label = "Analytical")
diff_prof = lines!(ax1, diff_1[:, 1], zw, color=:red3, linestyle=:dash, linewidth = 1.5,label = "M-A")
axislegend(ax1, position = :rb)
# Plot "Martin curve" for comparison

# POC_last = interior(POCt[end], 1, 1, :) 
# POC_flux = POC_last * (-wₛ*day)
# martin = POC_flux[grid.Nz]*((z[grid.Nz]+z₀)./(z.+z₀)).^martin_b

# lines!(ax1, martin, z, linewidth = 2,linestyle=:dash, label = "Martin curve")

############################################################
######################## Monthly ###########################
############################################################
ax2 = Axis(fig[1, 2], ylabel="z (m)", xlabel="POC conc (mol m⁻³)")
xlims!(ax2, -0.5, 2.5)
ylims!(ax2, -1000, 0)
lines!(ax2, POCn2, zw, linewidth = 3,label = "Modeled")
POC_prof2 = lines!(ax2, Ana_2[:, 1], zw, linewidth = 1.5,label = "Analytical")
diff_prof2 = lines!(ax2, diff_2[:, 1], zw, color=:red3, linestyle=:dash, linewidth = 1.5,label = "M-A")
axislegend(ax2, position = :rb)

record(fig, "POC_export_remin_1D.mp4", 366:nt, framerate=20) do nn
    n[] = nn
    POC_prof[1] = Ana_1[:, nn]
    POC_prof2[1] = Ana_2[:, nn]
    diff_prof[1] = diff_1[:, nn]
    diff_prof2[1] = diff_2[:, nn]
end
nothing