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

@inline κ(z, t) = 1e-4 + 1e-2 * exp(z / 25) + 1e-2 * exp(-(z + 1000) / 50)
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
                                                                                        # linear_remineralization_rate = 0.1/day,                                                    
                                                                                        maximum_bacteria_growth_rate = 0.9/day,
                                                                                        detritus_half_saturation = 0.05, 
                                                                                        detritus_vertical_velocity = -5/day,
                                                                                        stoichoimetric_ratio_iron_to_nitrate  = 2e-4),
                                    tracers = (:N, :P, :Z, :B, :D, :F),
                                    tracer_advection = WENO(),
                                    buoyancy = nothing,
                                    closure = vertical_diffusion) 

# ## Initial conditions
#
# We initialize the model with reasonable nutrients, detritus, and a nutrient
# mixed layer.

set!(model1, N=20, P=1e-1, Z=0, B=1e-1, D=5e-1, F=1e-3)

simulation1 = Simulation(model1, Δt=30minutes, stop_time=1000days)

# Define a callback function to apply the perturbation

# function perturbation_callback(sim)
#     # For a month in Year 5: 1810-1840 days
#     if time(sim) >= 2000days && time(sim) < 2030days 

#         # Apply Fe perturbation
#         Fₘ = model1.tracers.F 
#         Fₘ[:,:,95:100] .+=0.002
#         set!(model1, F=Fₘ)    
#     end
# end

# # Add the perturbation callback to the simulation
# simulation1.callbacks[:perturbation] = Callback(perturbation_callback,IterationInterval(100))

function progress(sim)
    @printf("Iteration: %d, time: %s, total(N): %.2e , total(Fe): %.2e\n",
            iteration(sim), prettytime(sim),
            sum(model1.tracers.N)+sum(model1.tracers.P)+sum(model1.tracers.Z)+sum(model1.tracers.B)+sum(model1.tracers.D),
            sum(model1.tracers.F))
    return nothing
end

simulation1.callbacks[:progress] = Callback(progress, IterationInterval(100))

N = model1.tracers.N
P = model1.tracers.P
Z = model1.tracers.Z
B = model1.tracers.B
D = model1.tracers.D
F = model1.tracers.F

filename = "NPBD_Fe.jld2"

simulation1.output_writers[:fields] = JLD2OutputWriter(model1, model1.tracers;
                                                      filename,
                                                      schedule = TimeInterval(1day),
                                                      overwrite_existing = true)

run!(simulation1)

#=
####################################### Part 2: Fe fertilization #######################################
Nₘ = model1.tracers.N
Pₘ = model1.tracers.P
Zₘ = model1.tracers.Z
Bₘ = model1.tracers.B
Dₘ = model1.tracers.D

model2 = HydrostaticFreeSurfaceModel(; grid,
                                    velocities = PrescribedVelocityFields(),
                                    biogeochemistry = NutrientsPlanktonBacteriaDetritus(; grid,
                                                                                        #linear_remineralization_rate = 0.1/day,                                                    
                                                                                        maximum_bacteria_growth_rate = 0.9/day,
                                                                                        detritus_half_saturation = 0.05, 
                                                                                        detritus_vertical_velocity = -5/day,
                                                                                        iron_concentration =0.0001),
                                    tracers = (:N, :P, :Z, :B, :D),
                                    tracer_advection = WENO(),
                                    buoyancy = nothing,
                                    closure = vertical_diffusion) 

set!(model2, N=Nₘ, P=Pₘ, Z=Zₘ,B=Bₘ, D=Dₘ)
simulation2 = Simulation(model2, Δt=30minutes, stop_time=100days)

function progress(sim)
    @printf("Iteration: %d, time: %s, total(N): %.2e \n",
            iteration(sim), prettytime(sim),
            sum(model2.tracers.N)+sum(model2.tracers.P)+sum(model2.tracers.Z)+sum(model2.tracers.B)+sum(model2.tracers.D))
    return nothing
end
simulation2.callbacks[:progress] = Callback(progress, IterationInterval(100))

filename = "NPBD_FeFert2.jld2"
simulation2.output_writers[:fields] = JLD2OutputWriter(model2, model2.tracers;
                                                      filename = filename,
                                                      schedule = TimeInterval(1day),
                                                      overwrite_existing = true)

run!(simulation2)

####################################### Part 3: After Fe fertilization #######################################
Nₘ2 = model2.tracers.N
Pₘ2 = model2.tracers.P
Zₘ2 = model2.tracers.Z
Bₘ2 = model2.tracers.B
Dₘ2 = model2.tracers.D

model3 = HydrostaticFreeSurfaceModel(; grid,
                                    velocities = PrescribedVelocityFields(),
                                    biogeochemistry = NutrientsPlanktonBacteriaDetritus(; grid,
                                                                                        # linear_remineralization_rate = 0.1/day,                                                    
                                                                                        maximum_bacteria_growth_rate = 0.9/day,
                                                                                        detritus_half_saturation = 0.05, 
                                                                                        detritus_vertical_velocity = -5/day,
                                                                                        iron_concentration =0.00003),
                                    tracers = (:N, :P, :Z, :B, :D),
                                    tracer_advection = WENO(),
                                    buoyancy = nothing,
                                    closure = vertical_diffusion) 

set!(model3, N=Nₘ2, P=Pₘ2, Z=Zₘ2,B=Bₘ2, D=Dₘ2)
simulation3 = Simulation(model3, Δt=30minutes, stop_time=1825days)

function progress(sim)
    @printf("Iteration: %d, time: %s, total(N): %.2e \n",
            iteration(sim), prettytime(sim),
            sum(model3.tracers.N)+sum(model3.tracers.P)+sum(model3.tracers.Z)+sum(model3.tracers.B)+sum(model3.tracers.D))
    return nothing
end
simulation3.callbacks[:progress] = Callback(progress, IterationInterval(100))

filename = "NPBD_FeFert3.jld2"
simulation3.output_writers[:fields] = JLD2OutputWriter(model3, model3.tracers;
                                                      filename = filename,
                                                      schedule = TimeInterval(1day),
                                                      overwrite_existing = true)

run!(simulation3)
=#

# ###################### Visualization ######################

# All that's left is to visualize the results.
Nt = FieldTimeSeries(filename, "N")
Pt = FieldTimeSeries(filename, "P")
Bt = FieldTimeSeries(filename, "B")
Dt = FieldTimeSeries(filename, "D")
Ft = FieldTimeSeries(filename, "F")

t = Pt.times
nt = length(t)
z = znodes(Pt)

fig = Figure(;size=(1200, 600))#resolution=(1200, 600))

axN = Axis(fig[1, 1], ylabel="z (m)", xlabel="[Nutrient] (mmol m⁻³)")
axB = Axis(fig[1, 2], ylabel="z (m)", xlabel="[Biomass] (mmol m⁻³)")
axD = Axis(fig[1, 3], ylabel="z (m)", xlabel="[Detritus] (mmol m⁻³)")
axF = Axis(fig[1, 4], ylabel="z (m)", xlabel="[Iron] (mmol m⁻³)")

xlims!(axN, -0.1, 50)
xlims!(axB, -0.1, 1.0)
xlims!(axD, -0.1, 1.5)
xlims!(axF, 0, 5e-3)

slider = Slider(fig[2, 1:4], range=1:nt, startvalue=1)
n = slider.value

title = @lift @sprintf("Equilibrium biogeochemistry at t = %d days", t[$n] / day)
Label(fig[0, 1:4], title)

Nn = @lift interior(Nt[$n], 1, 1, :)
Pn = @lift interior(Pt[$n], 1, 1, :)
Bn = @lift interior(Bt[$n], 1, 1, :)
Dn = @lift interior(Dt[$n], 1, 1, :)
Fn = @lift interior(Ft[$n], 1, 1, :)

lines!(axN, Nn, z, label="N")
lines!(axB, Pn, z, label="P")
lines!(axB, Bn, z, label="B")
lines!(axD, Dn, z, label="D")
lines!(axF, Fn, z, label="Fe")
axislegend(axB, position = :rb)

record(fig, "nutrients_plankton_bacteria_detritus.mp4", 1:nt, framerate=24) do nn
    n[] = nn
end
nothing
