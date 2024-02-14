# # Nutrients, plankton, bacteria, detritus
#
# This example illustrates how to use ClimaOceanBiogeochemistry's
# `NutrientsPlanktonBacteriaDetrius` model in a single column context.

using ClimaOceanBiogeochemistry: NutrientsPlanktonBacteriaDetritus

using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Printf
using CairoMakie

# ## A single column grid
#
# We set up a single column grid whose depth is `H` and with `Nz` points

H = 1000
z = (-H, 0)
Nz = 100

grid = RectilinearGrid(size = Nz; z, topology = (Flat, Flat, Bounded))

# A prescribed vertical tracer diffusivity
# 
# We define a tracer diffusivity that mixes a lot near the surface
# (in the top 50 m), and less down below.

@inline κ(x, y, z, t) = 1e-4 + 1e-2 * exp(z / 25) + 1e-2 * exp(-(z + 1000) / 50)
vertical_diffusion = VerticalScalarDiffusivity(; κ)

# We put the pieces together.
# The important line here is `biogeochemistry = NutrientsPlanktonBacteriaDetritus(; grid)`.
# We use all default parameters.

model = HydrostaticFreeSurfaceModel(; grid,
                                    velocities = PrescribedVelocityFields(),
                                    biogeochemistry = NutrientsPlanktonBacteriaDetritus(; grid),
                                    tracers = (:N, :P, :Z, :B, :D1, :D2),
                                    tracer_advection = WENO(),
                                    buoyancy = nothing,
                                    closure = vertical_diffusion)

# ## Initial conditions
#
# We initialize the model with reasonable nutrients, detritus, and a nutrient
# mixed layer.

set!(model, N=3, P=1e-1, Z=1e-1, B=1e-1, D1=8e-2, D2=2e-2)

simulation = Simulation(model, Δt=30minutes, stop_time=100days)

function progress(sim)
    @printf("Iteration: %d, time: %s, total(N): %.2e \n",
            iteration(sim), prettytime(sim),
            sum(model.tracers.N) + sum(model.tracers.P) + sum(model.tracers.B) + sum(model.tracers.D1) + sum(model.tracers.D2))
    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

# Let's see the initial condition

N = model.tracers.N
P = model.tracers.P
Z = model.tracers.Z
B = model.tracers.B
D1 = model.tracers.D1
D2 = model.tracers.D2

z = znodes(N)

fig = Figure()

axN  = Axis(fig[1, 1], xlabel="Nutrient concentration (N)", ylabel="z (m)")
axP  = Axis(fig[1, 2], xlabel="Phytoplankton concentration (P)", ylabel="z (m)")
axZ  = Axis(fig[1, 3], xlabel="Zooplankton concentration (Z)", ylabel="z (m)")
axB  = Axis(fig[1, 4], xlabel="Bacteria concentration (B)", ylabel="z (m)")
axD1 = Axis(fig[1, 5], xlabel="Detritus concentration (D1)", ylabel="z (m)")
axD2 = Axis(fig[1, 6], xlabel="Detritus concentration (D2)", ylabel="z (m)")

lines!(axN, interior(N, 1, 1, :), z)
lines!(axP, interior(P, 1, 1, :), z)
lines!(axZ, interior(Z, 1, 1, :), z)
lines!(axB, interior(B, 1, 1, :), z)
lines!(axD1, interior(D1, 1, 1, :), z)
lines!(axD2, interior(D2, 1, 1, :), z)

current_figure()

# Now we add an output writer to the simulation and run the simulation.

filename = "nutrients_plankton_bacteria_detritus.jld2"

simulation.output_writers[:fields] = JLD2OutputWriter(model, model.tracers;
                                                      filename,
                                                      schedule = TimeInterval(1day),
                                                      overwrite_existing = true)

run!(simulation)

# ## Visualization
#
# All that's left is to visualize the results.

Pt = FieldTimeSeries(filename, "P")
Zt = FieldTimeSeries(filename, "Z")
Bt = FieldTimeSeries(filename, "B")
D1t = FieldTimeSeries(filename, "D1")
D2t = FieldTimeSeries(filename, "D2")
Nt = FieldTimeSeries(filename, "N")

t = Pt.times
nt = length(t)
z = znodes(Pt)

fig = Figure(size=(1200, 600))

axN  = Axis(fig[1, 1], xlabel="[Nutrient] (mmol m⁻³)", ylabel="z (m)")
axP  = Axis(fig[1, 2], xlabel="[Phytoplankton] (mmol m⁻³)")
axZ  = Axis(fig[1, 3], xlabel="[Zooplankton] (mmol m⁻³)")
axB  = Axis(fig[1, 4], xlabel="[Bacteria] (mmol m⁻³)")
axD1 = Axis(fig[1, 5], xlabel="[Dissolved Detritus] (mmol m⁻³)")
axD2 = Axis(fig[1, 6], xlabel="[Particulate Detritus] (mmol m⁻³)")

slider = Slider(fig[2, 1:6], range=1:nt, startvalue=1)
n = slider.value

title = @lift @sprintf("Equilibrium biogeochemistry at t = %d days", t[$n] / day)
Label(fig[0, 1:6], title)

Nn  = @lift interior(Nt[$n], 1, 1, :)
Pn  = @lift interior(Pt[$n], 1, 1, :)
Zn  = @lift interior(Zt[$n], 1, 1, :)
Bn  = @lift interior(Bt[$n], 1, 1, :)
D1n = @lift interior(D1t[$n], 1, 1, :)
D2n = @lift interior(D2t[$n], 1, 1, :)

lines!(axP, Pn, z)
lines!(axZ, Zn, z)
lines!(axD1, D1n, z)
lines!(axD2, D2n, z)
lines!(axB, Bn, z)
lines!(axN, Nn, z)

record(fig, "nutrients_plankton_bacteria_detritus.mp4", 1:nt, framerate=24) do nn
    n[] = nn
end
nothing #hide

# ![](nutrients_plankton_bacteria_detritus.mp4)

# Let's plot a snapshot of the last frame.

Nn_last  = interior(Nt[end], 1, 1, :)
Pn_last  = interior(Pt[end], 1, 1, :)
Zn_last  = interior(Zt[end], 1, 1, :)
Bn_last  = interior(Bt[end], 1, 1, :)
D1n_last = interior(D1t[end], 1, 1, :)
D2n_last = interior(D2t[end], 1, 1, :)

last_frame = Figure(size=(1200, 600))
axN  = Axis(last_frame[1, 1], xlabel="[N] (mmol m⁻³)", ylabel="z (m)")
axP  = Axis(last_frame[1, 2], xlabel="[P] (mmol m⁻³)")
axZ  = Axis(last_frame[1, 3], xlabel="[Z] (mmol m⁻³)")
axB  = Axis(last_frame[1, 4], xlabel="[B] (mmol m⁻³)")
axD1 = Axis(last_frame[1, 5], xlabel="[dD] (mmol m⁻³)")
axD2 = Axis(last_frame[1, 6], xlabel="[pD] (mmol m⁻³)")

lines!(axP, Pn_last, z)
lines!(axZ, Zn_last, z)
lines!(axD1, D1n_last, z)
lines!(axD2, D2n_last, z)
lines!(axB, Bn_last, z)
lines!(axN, Nn_last, z)

# Save the last frame as a figure
save("NPZDB.png", last_frame)

# ![](NPZDB.png)

# Another figure: we plot the sum of each variable against time.

N_time  = zeros(1:nt)
P_time  = zeros(1:nt)
Z_time  = zeros(1:nt)
B_time  = zeros(1:nt)
D1_time = zeros(1:nt)
D2_time = zeros(1:nt)

for times = 1:nt
    N_time[times] = sum(Nt[:, :, :, times])
    P_time[times] = sum(Pt[:, :, :, times])
    Z_time[times] = sum(Zt[:, :, :, times])
    B_time[times] = sum(Bt[:, :, :, times])
    D1_time[times] = sum(D1t[:, :, :, times])
    D2_time[times] = sum(D2t[:, :, :, times])
end

TimeVar = Figure()
ax2 = Axis(TimeVar[1, 1], title="Variable over time", ylabel="Variable (mmol m⁻³)", xlabel="Time (days)")
lines!(ax2, 1:nt, N_time, label="N")  
lines!(ax2, 1:nt, P_time, label="P")  
lines!(ax2, 1:nt, Z_time, label="Z")  
lines!(ax2, 1:nt, B_time, label="B")  
lines!(ax2, 1:nt, D1_time, label="dD") 
lines!(ax2, 1:nt, D2_time, label="pD") 

axislegend()

save("TimeVariations.png", TimeVar)

# ![](TimeVariations.png)
