# Build a one-dimensional single column model, 
# simulating POC remineralization and sinking in a minimalistic manner
# Purpose: Test the remineralization formulation of POC as a function of depth

using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Fields: ConstantField, ZeroField
using Oceananigans.Grids: Center, znode
using CairoMakie
using Printf

######################################### Grid ##########################################
# the 1D domain
H = 1000 # deepest depth 
Nz = 1000 # number of vertical grids (resolution = H/Nz)

grid = RectilinearGrid(size = Nz; z = (-H, 0), topology = (Flat, Flat, Bounded))

################################### Boundary condition ###################################
# Add a top boundary condition (fixed value) to ensure there is POC production from the top 
top_value = 5
value_top_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(top_value))

####################################### Forcing #######################################
# Particle forcing: remineralization and sinking

# Calculate remineralization of particulate organic phosphorus according to a first-order rate constant.
z₀ = -5 # Reference depth (m)
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
set!(model, POC=5)

simulation = Simulation(model, Δt=10minutes, stop_time=120days)

function progress(sim)
    @printf("Iteration: %d, time: %s, total(POC): %.2e \n",
            iteration(sim), prettytime(sim),
            sum(model.tracers.POC))
    return nothing
end
simulation.callbacks[:progress] = Callback(progress, IterationInterval(1000))

############################### Output file ###############################
filename = "POC_minimal.jld2"
simulation.output_writers[:fields] = JLD2OutputWriter(model, model.tracers;
                                                      filename,
                                                      schedule = TimeInterval(1day),
                                                      overwrite_existing = true)

# Start simulation!
run!(simulation)

# ###################### Visualization ######################

# All that's left is to visualize the results.
POCt = FieldTimeSeries(filename, "POC")

t = POCt.times
nt = length(t)
z = znodes(POCt)

# 1. Movie
fig = Figure(;size=(800, 600))

ax1 = Axis(fig[1, 1], ylabel="z (m)", xlabel="POC flux (mol m⁻² d⁻¹)")
ax2 = Axis(fig[1, 2], ylabel="z (m)", xlabel="Remin rate (mol m⁻³ d⁻¹)")

xlims!(ax1, -5, 55)
ylims!(ax1, -1000, 0)

xlims!(ax2, -0.5, 5)
ylims!(ax2, -1000, 0)

slider = Slider(fig[2, 1:2], range=1:nt, startvalue=1)
n = slider.value
title = @lift @sprintf("t = %d days", t[$n] / day)
Label(fig[0, 1:2], title)

# tracer concentration
POCn = @lift interior(POCt[$n], 1, 1, :)
# particle flux = concentration x sinking speed
FPOC = @lift (-wₛ*day)*(interior(POCt[$n], 1, 1, :))

lines!(ax1, FPOC, z, linewidth = 1.5,label = "Modeled flux")

# Plot "Martin curve" for comparison

POC_last = interior(POCt[end], 1, 1, :) 
POC_flux = POC_last * (-wₛ*day)
martin = POC_flux[grid.Nz]*((z[grid.Nz]+z₀)./(z.+z₀)).^martin_b

lines!(ax1, martin, z, linewidth = 2,linestyle=:dash, label = "Martin curve")
axislegend(ax1, position = :rb)

# 2. Plot rates
ReminRate = @lift (martin_b*(-wₛ*day)./(-z.+z₀)).*(interior(POCt[$n], 1, 1, :))
lines!(ax2, ReminRate, z,linewidth = 1.5,label = "Remin rate")

record(fig, "POC_minimal.mp4", 1:nt, framerate=24) do nn
    n[] = nn
end
nothing