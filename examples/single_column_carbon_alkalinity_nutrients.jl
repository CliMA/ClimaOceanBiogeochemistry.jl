# # Nutrients, plankton, bacteria, detritus
#
# This example illustrates how to use ClimaOceanBiogeochemistry's
# `NutrientsPlanktonBacteriaDetrius` model in a single column context.

using ClimaOceanBiogeochemistry: NutrientsPlanktonBacteriaDetritus, CarbonAlkalinityNutrients

using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Printf
using CairoMakie

# ## A single column grid
#
# We set up a single column grid with 4 m grid spacing that's 256 m deep:

grid = RectilinearGrid(size = 64,
                       z = (-256meters, 0),
                       topology = (Flat, Flat, Bounded))

# ## Convection then quiet
#
# To illustrate the dynamics of `NutrientsPlanktonBacteriaDetritus`,
# we set up a physical scenario in which strong convection drives turbulent mixing
# for 4 days, and then abruptly shuts off. Once the convective turbulence dies
# down, plankton start to grow.

Qᵇ(x, y, t) = ifelse(t < 4days, 1e-7, 0.0)
b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵇ))

# We put the pieces together.
# The important line here is `biogeochemistry = NutrientsPlanktonBacteriaDetritus()`.
# We use all default parameters.

model = HydrostaticFreeSurfaceModel(; grid,
                                    biogeochemistry = CarbonAlkalinityNutrients(),
                                    closure = CATKEVerticalDiffusivity(),
                                    tracers = (:b, :e),
                                    buoyancy = BuoyancyTracer(),
                                    boundary_conditions = (; b=b_bcs))

# ## Initial conditions
#
# We initialize the model with reasonable nutrients, detritus, and a nutrient
# mixed layer.

N₀ = 1e-1 # Surface nutrient concentration
D₀ = 1e-1 # Surface detritus concentration
dᴺ = 50.0 # Nutrient mixed layer depth
N² = 1e-5 # Buoyancy gradient, s⁻²

bᵢ(x, y, z) = N² * z
Nᵢ(x, y, z) = N₀ * max(1, exp(-(z + dᴺ) / 100))
Dᵢ(x, y, z) = D₀ * exp(z / 10)

#set!(model, b=bᵢ, P=1e-1, B=1e-1, D=Dᵢ, N=Nᵢ, e=1e-6)
set!(model, b=bᵢ, e=1e-6)

# ## A simulation of physical-biological interaction
# 
# We construct a simple simulation that emits a message every 10 iterations
# and outputs tracer fields.

simulation = Simulation(model, Δt=10minutes, stop_time=12days)

progress(sim) = @printf("Iteration: %d, time: %s\n", iteration(sim), prettytime(sim))
simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

filename = "single_column_carbon_alkalinity_nutrients.jld2"

simulation.output_writers[:fields] = JLD2OutputWriter(model, model.tracers;
                                                      filename,
                                                      schedule = TimeInterval(20minutes),
                                                      overwrite_existing = true)

run!(simulation)

# ## Visualization
#
# All that's left is to visualize the results.

bt   = FieldTimeSeries(filename, "b")
et   = FieldTimeSeries(filename, "e")
DICt = FieldTimeSeries(filename, "DIC")
Alkt = FieldTimeSeries(filename, "Alk")
PO₄t = FieldTimeSeries(filename, "PO₄")
NO₃t = FieldTimeSeries(filename, "NO₃")
DOPt = FieldTimeSeries(filename, "DOP")
Fet  = FieldTimeSeries(filename, "Fe")

t = bt.times
nt = length(t)
z = znodes(bt)

fig = Figure(resolution=(1200, 600))

axb = Axis(fig[1, 1], ylabel="z (m)", xlabel="Buoyancy (m² s⁻³)")
axe = Axis(fig[1, 2], ylabel="z (m)", xlabel="Turbulent kinetic energy (m² s²)")
axP = Axis(fig[1, 3], ylabel="z (m)", xlabel="Concentration (mmol)")
axN = Axis(fig[1, 4], ylabel="z (m)", xlabel="Nutrient concentration (mmol)")

xlims!(axe, -1e-5, 1e-3)
xlims!(axP, 0, 0.2)

slider = Slider(fig[2, 1:4], range=1:nt, startvalue=1)
n = slider.value

title = @lift @sprintf("Convecting plankton at t = %d hours", t[$n] / hour)
Label(fig[0, 1:4], title)

bn = @lift interior(bt[$n], 1, 1, :)
en = @lift interior(et[$n], 1, 1, :)
DICn = @lift interior(DICt[$n], 1, 1, :) 
Alkn = @lift interior(Alkt[$n], 1, 1, :) 
PO₄n = @lift interior(PO₄t[$n], 1, 1, :) 
NO₃n = @lift interior(NO₃t[$n], 1, 1, :) 
DOPn = @lift interior(DOPt[$n], 1, 1, :) 
Fen  = @lift interior(Fet[$n], 1, 1, :)  

lines!(axb, bn, z)
lines!(axe, en, z)

lines!(axP, DICn, z, label="DIC")
lines!(axP, PO₄n, z, label="Phosphate")
lines!(axP, Fen, z, label="Iron")
axislegend(axP)

lines!(axN, DOPn, z)

record(fig, "carbon_alkalinity_nutrients.mp4", 1:nt, framerate=24) do nn
    n[] = nn
end
nothing #hide

# ![](carbon_alkalinity_nutrients.mp4)
#
