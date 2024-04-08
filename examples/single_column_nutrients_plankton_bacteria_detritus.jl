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

# We put the pieces together.
# The important line here is `biogeochemistry = NutrientsPlanktonBacteriaDetritus()`.
# We use all default parameters.

model = HydrostaticFreeSurfaceModel(; grid,
                                    velocities = PrescribedVelocityFields(),
                                    biogeochemistry = NutrientsPlanktonBacteriaDetritus(grid;
                                                                                        maximum_bacteria_growth_rate = 2/day,
                                                                                        detritus_half_saturation = 1.0),
                                                                                        #detritus_vertical_velocity = 0/day),
                                    tracers = (:N, :P, :Z, :B, :D),
                                    tracer_advection = WENO(),
                                    buoyancy = nothing,
                                    closure = vertical_diffusion) 
                                    #closure = CATKEVerticalDiffusivity(),
                                    #tracers = (:b, :e),
                                    #buoyancy = BuoyancyTracer(),
                                    #boundary_conditions = (; b=b_bcs))

# ## Initial conditions
#
# We initialize the model with reasonable nutrients, detritus, and a nutrient
# mixed layer.

#set!(model, N1=3, P1=1e-1, P2=5e-2,Z1=1e-1,Z2=1e-2, B1=1e-1,B2=2e-2, D1=1e-1, D2=1e-2,D3=2e-1, D4=1e-1,D5=2e-1)
set!(model, N=20, P=1e-1, Z=0,B=1e-1, D=1e-1)

simulation = Simulation(model, Δt=30minutes, stop_time=3650days)

function progress(sim)
    @printf("Iteration: %d, time: %s, total(N): %.3e \n",
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

simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))

N = model.tracers.N
P = model.tracers.P
Z = model.tracers.Z
B = model.tracers.B
D = model.tracers.D

#=
N1 = model.tracers.N1
P1 = model.tracers.P1
P2 = model.tracers.P2
Z1 = model.tracers.Z1
Z2 = model.tracers.Z2
B1 = model.tracers.B1
B2 = model.tracers.B2
D1 = model.tracers.D1 
D2 = model.tracers.D2
D3 = model.tracers.D3
D4 = model.tracers.D4
D5 = model.tracers.D5
=#

z = znodes(N)

fig = Figure()

axN = Axis(fig[1, 1], xlabel="Nutrient concentration (N)", ylabel="z (m)")
axP = Axis(fig[1, 2], xlabel="Phytoplankton concentration (P)", ylabel="z (m)")
#axZ = Axis(fig[1, 3], xlabel="Zooplankton concentration (Z)", ylabel="z (m)")
axB = Axis(fig[1, 3], xlabel="Bacteria concentration (B)", ylabel="z (m)")
axD = Axis(fig[1, 4], xlabel="Detritus concentration (D)", ylabel="z (m)")

#=
lines!(axN, interior(N1, 1, 1, :), z)
lines!(axP, interior(P1, 1, 1, :), z)
lines!(axP, interior(P2, 1, 1, :), z)
lines!(axZ, interior(Z1, 1, 1, :), z)
lines!(axZ, interior(Z2, 1, 1, :), z)
lines!(axB, interior(B1, 1, 1, :), z)
lines!(axB, interior(B2, 1, 1, :), z)
lines!(axD, interior(D1, 1, 1, :), z)
lines!(axD, interior(D2, 1, 1, :), z)
lines!(axD, interior(D3, 1, 1, :), z)
lines!(axD, interior(D4, 1, 1, :), z)
lines!(axD, interior(D5, 1, 1, :), z)
=#
lines!(axN, interior(N, 1, 1, :), z)
lines!(axP, interior(P, 1, 1, :), z)
lines!(axB, interior(B, 1, 1, :), z)
lines!(axD, interior(D, 1, 1, :), z)

display(fig)

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

#=
P1t = FieldTimeSeries(filename, "P1")
Z1t = FieldTimeSeries(filename, "Z1")
B1t = FieldTimeSeries(filename, "B1")
P2t = FieldTimeSeries(filename, "P2")
Z2t = FieldTimeSeries(filename, "Z2")
B2t = FieldTimeSeries(filename, "B2")

D1t = FieldTimeSeries(filename, "D1")
D2t = FieldTimeSeries(filename, "D2")
D3t = FieldTimeSeries(filename, "D3")
D4t = FieldTimeSeries(filename, "D4")
D5t = FieldTimeSeries(filename, "D5")
N1t = FieldTimeSeries(filename, "N1")
=#

t = Pt.times
nt = length(t)
z = znodes(Pt)

fig = Figure(;size=(1200, 600))#resolution=(1200, 600))

axN = Axis(fig[1, 1], ylabel="z (m)", xlabel="[Nutrient] (mmol m⁻³)")
axP = Axis(fig[1, 2], ylabel="z (m)", xlabel="[Phytoplankton] (mmol m⁻³)")
#axZ = Axis(fig[1, 3], ylabel="z (m)", xlabel="[Zooplankton] (mmol m⁻³)")
axB = Axis(fig[1, 3], ylabel="z (m)", xlabel="[Bacteria] (mmol m⁻³)")
axD = Axis(fig[1, 4], ylabel="z (m)", xlabel="[Detritus] (mmol m⁻³)")

xlims!(axN, -0.1, 30)
xlims!(axP, -0.1, 1.5)
xlims!(axB, -0.1, 1.5)
xlims!(axD, -0.1, 20)

slider = Slider(fig[2, 1:4], range=1:nt, startvalue=1)
n = slider.value

title = @lift @sprintf("Equilibrium biogeochemistry at t = %d days", t[$n] / day)
Label(fig[0, 1:4], title)

#=
N1n = @lift interior(N1t[$n], 1, 1, :)
P1n = @lift interior(P1t[$n], 1, 1, :)
P2n = @lift interior(P2t[$n], 1, 1, :)
Z1n = @lift interior(Z1t[$n], 1, 1, :)
Z2n = @lift interior(Z2t[$n], 1, 1, :)
B1n = @lift interior(B1t[$n], 1, 1, :)
B2n = @lift interior(B2t[$n], 1, 1, :)
D1n = @lift interior(D1t[$n], 1, 1, :)
D2n = @lift interior(D2t[$n], 1, 1, :)
D3n = @lift interior(D3t[$n], 1, 1, :)
D4n = @lift interior(D4t[$n], 1, 1, :)
D5n = @lift interior(D5t[$n], 1, 1, :)

lines!(axN, N1n, z)
lines!(axP, P1n, z)
lines!(axZ, Z1n, z)
lines!(axP, P2n, z)
lines!(axZ, Z2n, z)
lines!(axB, B1n, z)
lines!(axB, B2n, z)
lines!(axD, D1n, z)
lines!(axD, D2n, z)
=#

Nn = @lift interior(Nt[$n], 1, 1, :)
Pn = @lift interior(Pt[$n], 1, 1, :)
Bn = @lift interior(Bt[$n], 1, 1, :)
Dn = @lift interior(Dt[$n], 1, 1, :)

lines!(axN, Nn, z)
lines!(axP, Pn, z)
lines!(axB, Bn, z)
lines!(axD, Dn, z)

record(fig, "nutrients_plankton_bacteria_detritus.mp4", 1:nt, framerate=24) do nn
    n[] = nn
end
nothing #hide

# ![](nutrients_plankton_bacteria_detritus.mp4)

#= Figure 1: Extract the last frame
N1n_last = interior(N1t[end], 1, 1, :)
P1n_last = interior(P1t[end], 1, 1, :)
P2n_last = interior(P2t[end], 1, 1, :)
Z1n_last = interior(Z1t[end], 1, 1, :)
Z2n_last = interior(Z2t[end], 1, 1, :)
B1n_last = interior(B1t[end], 1, 1, :)
B2n_last = interior(B2t[end], 1, 1, :)
D1n_last = interior(D1t[end], 1, 1, :)
D2n_last = interior(D2t[end], 1, 1, :)
D3n_last = interior(D3t[end], 1, 1, :)
D4n_last = interior(D4t[end], 1, 1, :)
D5n_last = interior(D5t[end], 1, 1, :)

last_frame = Figure(resolution=(1200, 600))
axN = Axis(last_frame[1, 1], ylabel="z (m)", xlabel="[N] (mmol m⁻³)")
axP = Axis(last_frame[1, 2], xlabel="[P] (mmol m⁻³)")
axZ = Axis(last_frame[1, 3], xlabel="[Z] (mmol m⁻³)")
axB = Axis(last_frame[1, 4], xlabel="[B] (mmol m⁻³)")
axD = Axis(last_frame[1, 5], xlabel="[D] (mmol m⁻³)")

lines!(axP, P1n_last, z)
lines!(axP, P2n_last, z)
lines!(axZ, Z1n_last, z)
lines!(axZ, Z2n_last, z)
lines!(axD, D1n_last, z)
lines!(axD, D2n_last, z)
lines!(axD, D3n_last, z)
lines!(axD, D4n_last, z)
lines!(axD, D5n_last, z)
lines!(axB, B1n_last, z)
lines!(axB, B2n_last, z)
lines!(axN, N1n_last, z)

# Save the last frame as a figure
save("test.png", last_frame)

#=
# Figure 2: sum of each variable vs. time
N_time = zeros(1:nt)
P_time = zeros(1:nt)
Z_time = zeros(1:nt)
B_time = zeros(1:nt)
D1_time = zeros(1:nt)
D2_time = zeros(1:nt)

for times = 1:nt
    N_time[times] = sum(Nt[:,:,:,times])
    P_time[times] = sum(Pt[:,:,:,times])
    Z_time[times] = sum(Zt[:,:,:,times])
    B_time[times] = sum(Bt[:,:,:,times])
    D1_time[times] = sum(D1t[:,:,:,times])
    D2_time[times] = sum(D2t[:,:,:,times])
end

TimeVar = Figure()
ax2 = Axis(TimeVar[1,1], title="Variable over time",ylabel="Variable (mmol m⁻³)", xlabel="Time (days)")
lines!(ax2, 1:nt, N_time, label="N")  
lines!(ax2, 1:nt, P_time, label="P")  
lines!(ax2, 1:nt, Z_time, label="Z")  
lines!(ax2, 1:nt, B_time, label="B")  
lines!(ax2, 1:nt, D1_time, label="dD") 
lines!(ax2, 1:nt, D2_time, label="pD") 

axislegend()
 
save("TimeVariations.png", TimeVar)
=#
=#