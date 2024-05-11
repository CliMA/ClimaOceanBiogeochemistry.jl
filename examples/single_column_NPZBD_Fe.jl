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
 
model = HydrostaticFreeSurfaceModel(; grid,
                                    velocities = PrescribedVelocityFields(),
                                    biogeochemistry = NutrientsPlanktonBacteriaDetritus(; grid,
                                                                                        #linear_remineralization_rate = 0.1/day,                                                    
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

#set!(model, N1=3, P1=1e-1, P2=5e-2,Z1=1e-1,Z2=1e-2, B1=1e-1,B2=2e-2, D1=1e-1, D2=1e-2,D3=2e-1, D4=1e-1,D5=2e-1)
set!(model, N=20, P=1e-1, Z=0,B=1e-1, D=5e-1)

simulation = Simulation(model, Δt=30minutes, stop_time=2000days)

function progress(sim)
    @printf("Iteration: %d, time: %s, total(N): %.2e \n",
            iteration(sim), prettytime(sim),
            sum(model.tracers.N)+sum(model.tracers.P)+sum(model.tracers.Z)+sum(model.tracers.B)+sum(model.tracers.D))
            #=
            sum(model.tracers.N1)+
            sum(model.tracers.P1)+sum(model.tracers.P2)+
            sum(model.tracers.Z1)+ sum(model.tracers.Z2)
            sum(model.tracers.B1)+ sum(model.tracers.B2)+
            sum(model.tracers.D1)+sum(model.tracers.D2)+sum(model.tracers.D3)+sum(model.tracers.D4)+sum(model.tracers.D5))
            =#
    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

N = model.tracers.N
P = model.tracers.P
Z = model.tracers.Z
B = model.tracers.B
D = model.tracers.D

filename = "nutrients_plankton_bacteria_detritus.jld2"

simulation.output_writers[:fields] = JLD2OutputWriter(model, model.tracers;
                                                      filename,
                                                      schedule = TimeInterval(1day),
                                                      overwrite_existing = true)

run!(simulation)

# ###################### Visualization ######################
#
# All that's left is to visualize the results.
Nt = FieldTimeSeries(filename, "N")
Pt = FieldTimeSeries(filename, "P")
Bt = FieldTimeSeries(filename, "B")
Dt = FieldTimeSeries(filename, "D")

t = Pt.times
nt = length(t)
z = znodes(Pt)

fig = Figure(;size=(1200, 600))#resolution=(1200, 600))

axN = Axis(fig[1, 1], ylabel="z (m)", xlabel="[Nutrient] (mmol m⁻³)")
axB = Axis(fig[1, 2], ylabel="z (m)", xlabel="[Biomass] (mmol m⁻³)")
axD = Axis(fig[1, 3], ylabel="z (m)", xlabel="[Detritus] (mmol m⁻³)")

xlims!(axN, -0.1, 50)
xlims!(axB, -0.1, 1.0)
xlims!(axD, -0.1, 1.5)

slider = Slider(fig[2, 1:3], range=1:nt, startvalue=1)
n = slider.value

title = @lift @sprintf("Equilibrium biogeochemistry at t = %d days", t[$n] / day)
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