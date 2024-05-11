# # Nutrients, plankton, bacteria, detritus
#
# This example illustrates how to use ClimaOceanBiogeochemistry's
# `NutrientsPlanktonBacteriaDetrius` model in a single column context.

using ClimaOceanBiogeochemistry: NutrientsPlanktonBacteriaDetritus

using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Printf
using GLMakie
#using CairoMakie

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

@inline κ(z, t) = 1e-4 + 1e-2 * exp(z / 50) + 1e-2 * exp(-(z + 1000) / 50)
vertical_diffusion = VerticalScalarDiffusivity(; κ)

# The following two lines should be later debugged in Oceananigans
import Oceananigans.Models.HydrostaticFreeSurfaceModels: validate_tracer_advection, AbstractAdvectionScheme, SingleColumnGrid
validate_tracer_advection(tracer_advection::AbstractAdvectionScheme, ::SingleColumnGrid) = tracer_advection, NamedTuple()

# We put the pieces together.
# The important line here is `biogeochemistry = NutrientsPlanktonBacteriaDetritus()`.
# We use all default parameters.
 
model1 = HydrostaticFreeSurfaceModel(; grid,
                                    velocities = PrescribedVelocityFields(),
                                    biogeochemistry = NutrientsPlanktonBacteriaDetritus(; grid,
                                                                                        #linear_remineralization_rate = 0.08/day,                                                    
                                                                                        maximum_bacteria_growth_rate = 0.9/day,
                                                                                        detritus_half_saturation = 0.05, 
                                                                                        detritus_vertical_velocity = -5/day),
                                    tracers = (:N, :P, :Z, :B, :D),
                                    tracer_advection = WENO(),
                                    buoyancy = nothing,
                                    closure = vertical_diffusion) 

# ## Initial conditions
#
# We initialize the model with reasonable nutrients, detritus, and a nutrient
# mixed layer.

set!(model1, N=20, P=1e-1, Z=0,B=1e-1, D=5e-1)

simulation1 = Simulation(model1, Δt=30minutes, stop_time=1825days)

# To illustrate the dynamics of `NutrientsPlanktonBacteriaDetritus`,
# we set up a scenario: Model reaches quasi-equilibrium state after 5 years (around 1826 days),
# then we add a perturbation to the system (e.g. a surface nutrient source)
# Check how dynamics would change in this ecosystem.

# Define a callback function to apply the perturbation

#=
function perturbation_callback(sim)
    # For a month in Year 5: 1810-1840 days
    if time(sim) >= 100days && time(sim) < 150days 
        
        # Apply MLD perturbation
        @inline κ(z, t) = 1e-4 + 1e-2 * exp(z / 25) + 1e-2 * exp(-(z + 1000) / 50)
        vertical_diffusion = VerticalScalarDiffusivity(; κ)

        set!(model, model.closure = vertical_diffusion)      
    end
end

# Add the perturbation callback to the simulation
simulation.callbacks[:perturbation] = Callback(perturbation_callback,IterationInterval(100))
=#

function progress(sim)
    @printf("Iteration: %d, time: %s, total(N): %.2e \n",
            iteration(sim), prettytime(sim),
            sum(model1.tracers.N)+sum(model1.tracers.P)+sum(model1.tracers.Z)+sum(model1.tracers.B)+sum(model1.tracers.D))
    return nothing
end

simulation1.callbacks[:progress] = Callback(progress, IterationInterval(100))

simulation1.output_writers[:fields] = JLD2OutputWriter(model1, model1.tracers;
                                                      filename = "NPBD_MLD2_before.jld2",
                                                      schedule = TimeInterval(1day),
                                                      overwrite_existing = true)

run!(simulation1)

####################################### Part 2: MLD shoaling #######################################
Nₘ = model1.tracers.N
Pₘ = model1.tracers.P
Zₘ = model1.tracers.Z
Bₘ = model1.tracers.B
Dₘ = model1.tracers.D

@inline κ(z, t) = 1e-4 + 1e-2 * exp(z / 75) + 1e-2 * exp(-(z + 1000) / 50)
vertical_diffusion = VerticalScalarDiffusivity(; κ)
 
model = HydrostaticFreeSurfaceModel(; grid,
                                    velocities = PrescribedVelocityFields(),
                                    biogeochemistry = NutrientsPlanktonBacteriaDetritus(; grid,
                                                                                        #linear_remineralization_rate = 0.08/day,                                                    
                                                                                        maximum_bacteria_growth_rate = 0.9/day,
                                                                                        detritus_half_saturation = 0.05, 
                                                                                        detritus_vertical_velocity = -5/day),
                                    tracers = (:N, :P, :Z, :B, :D),
                                    tracer_advection = WENO(),
                                    buoyancy = nothing,
                                    closure = vertical_diffusion) 

set!(model, N=Nₘ, P=Pₘ, Z=Zₘ,B=Bₘ, D=Dₘ)
simulation = Simulation(model, Δt=30minutes, stop_time=1825days)

function progress(sim)
    @printf("Iteration: %d, time: %s, total(N): %.2e \n",
            iteration(sim), prettytime(sim),
            sum(model.tracers.N)+sum(model.tracers.P)+sum(model.tracers.Z)+sum(model.tracers.B)+sum(model.tracers.D))
    return nothing
end
simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

filename = "NPBD_MLD2_after.jld2"
simulation.output_writers[:fields] = JLD2OutputWriter(model, model.tracers;
                                                      filename,
                                                      schedule = TimeInterval(1day),
                                                      overwrite_existing = true)

run!(simulation)

# ###################### Visualization ######################

Nt = FieldTimeSeries(filename, "N")
Pt = FieldTimeSeries(filename, "P")
Bt = FieldTimeSeries(filename, "B")
Dt = FieldTimeSeries(filename, "D")

t = Pt.times
nt = length(t)
z = znodes(Pt)

fig = Figure(;size=(1200, 600))

axN = Axis(fig[1, 1], ylabel="z (m)", xlabel="[Nutrient] (mmol m⁻³)")
axB = Axis(fig[1, 2], ylabel="z (m)", xlabel="[Biomass] (mmol m⁻³)")
axD = Axis(fig[1, 3], ylabel="z (m)", xlabel="[Detritus] (mmol m⁻³)")

xlims!(axN, -0.1, 50)
xlims!(axB, -0.1, 1.5)
xlims!(axD, -0.1, 3.0)

slider = Slider(fig[2, 1:3], range=1:(n0t+nt), startvalue=1)
n = slider.value

title = @lift @sprintf("Equilibrium biogeochemistry at t = %d days", t[$n] / day + 1825)
Label(fig[0, 1:3], title)

Nn = @lift interior(Nt[$n], 1, 1, :)
Pn = @lift interior(Pt[$n], 1, 1, :)
Bn = @lift interior(Bt[$n], 1, 1, :)
Dn = @lift interior(Dt[$n], 1, 1, :)

lines!(axN, Nn, z)
lines!(axB, Pn, z, label="P")
lines!(axB, Bn, z, label="B")
lines!(axD, Dn, z)

axislegend(axB, position = :rb)

record(fig, "nutrients_plankton_bacteria_detritus.mp4", 1:nt, framerate=24) do nn
    n[] = nn
end
