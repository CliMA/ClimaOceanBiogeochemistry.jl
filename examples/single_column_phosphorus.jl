# The phosphorus cycle: PO₄, DOP, POP
#
# This example illustrates how to use ClimaOceanBiogeochemistry's
# `CarbonAlkalinityNutrients` model in a single column context.

using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients

using Oceananigans
using Oceananigans.Units
using Printf
using CairoMakie

# ## A single column grid
#
# We set up a single column grid with Nz m grid spacing that's H m deep:

H = 1000
z = (-H, 0)
Nz = 100

grid = RectilinearGrid(size = Nz; z, topology = (Flat, Flat, Bounded))

# A prescribed vertical tracer diffusivity
# 
# We define a tracer diffusivity that mixes a lot near the surface
# (in the top 50 m), and less down below.

@inline κ(z, t) = 1e-4 + 1e-2 * exp(z / 25) + 1e-2 * exp(-(z + 1000) / 50)
vertical_diffusion = VerticalScalarDiffusivity(; κ)

# The following two lines should be later debugged in Oceananigans
import Oceananigans.Models.HydrostaticFreeSurfaceModels: validate_tracer_advection, AbstractAdvectionScheme, SingleColumnGrid
validate_tracer_advection(tracer_advection::AbstractAdvectionScheme, ::SingleColumnGrid) = tracer_advection, NamedTuple()

# We put the pieces together.
# The important line here is `biogeochemistry = CarbonAlkalinityNutrients()`.
# We use all default parameters.
 
model = HydrostaticFreeSurfaceModel(; grid,
                                    velocities = PrescribedVelocityFields(),
                                    biogeochemistry = CarbonAlkalinityNutrients(; grid,
                                                                                maximum_net_community_production_rate  = 2e-1/day,
                                                                                dissolved_organic_phosphorus_remin_timescale  = 1/365.25days,
                                                                                particulate_organic_phosphorus_sinking_velocity = -1/day),
                                    tracers = (:DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe),
                                    tracer_advection = WENO(),
                                    buoyancy = nothing,
                                    closure = vertical_diffusion) 

# ## Initial conditions
#
# We initialize the model with reasonable concentrations

set!(model, DIC=2.1, ALK=2.35, NO₃=3.6e-2, PO₄=2e-3, DOP=5e-4, POP=1.5e-3, Fe = 1e-6) # mol PO₄ m⁻³

simulation = Simulation(model, Δt=30minutes, stop_time=1000days)

function progress(sim)
    @printf("Iteration: %d, time: %s, total(P): %.2e \n",
            iteration(sim), prettytime(sim),
            sum(model.tracers.PO₄)+sum(model.tracers.DOP)+sum(model.tracers.POP))
    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

filename = "po4_dop_pop.jld2"

simulation.output_writers[:fields] = JLD2OutputWriter(model, model.tracers;
                                                      filename,
                                                      schedule = TimeInterval(1day),
                                                      overwrite_existing = true)

run!(simulation)

# ###################### Visualization ######################

# All that's left is to visualize the results.
PO4t = FieldTimeSeries(filename, "PO₄")
DOPt = FieldTimeSeries(filename, "DOP")
POPt = FieldTimeSeries(filename, "POP")

t = PO4t.times
nt = length(t)
z = znodes(PO4t)

fig = Figure(;size=(800, 600))#resolution=(1200, 600))

ax1 = Axis(fig[1, 1], ylabel="z (m)", xlabel="[Inorganic PO₄] (mol m⁻³)")
ax2 = Axis(fig[1, 2], ylabel="z (m)", xlabel="[Organic P] (mol m⁻³)")

xlims!(ax1, 0, 0.01)
xlims!(ax2, 0, 0.005)

slider = Slider(fig[2, 1:2], range=1:nt, startvalue=1)
n = slider.value

title = @lift @sprintf("Equilibrium biogeochemistry at t = %d days", t[$n] / day)
Label(fig[0, 1:2], title)

PO4n = @lift interior(PO4t[$n], 1, 1, :)
DOPn = @lift interior(DOPt[$n], 1, 1, :)
POPn = @lift interior(POPt[$n], 1, 1, :)

# martin = [NaN for i in 1:length(z)]
# for i in 1:length(z)
#     if z[i]<-115
#         martin[i] = exp.(z[i]/(-115)).^-0.84
#     end
# end

lines!(ax1, PO4n, z)
lines!(ax2, DOPn, z, label = "DOP")
lines!(ax2, POPn, z, label = "POP")
axislegend(ax2, position = :rb)

record(fig, "po4_dop_pop.mp4", 1:nt, framerate=24) do nn
    n[] = nn
end
nothing
