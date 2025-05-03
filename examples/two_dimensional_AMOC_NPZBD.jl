using GLMakie
# using CUDA
using Printf
using Statistics

using Oceananigans
using Oceananigans.Units
using Oceananigans.Fields: ZeroField, CenterField
using Oceananigans.BoundaryConditions: fill_halo_regions!

using ClimaOceanBiogeochemistry: NutrientsPlanktonBacteriaDetritus
using Oceananigans.Models.HydrostaticFreeSurfaceModels:
                    HydrostaticFreeSurfaceModel,
                    PrescribedVelocityFields
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization
using Oceananigans: TendencyCallsite

Ny = 500 
Nz = 200
Ly = 15000kilometers   # m
Lz = 4000           # m

arch = CPU()
# We use a two-dimensional grid, with a `Flat` `y`-direction:
grid = RectilinearGrid(arch,
                       size = (Ny, Nz),
                       y = (0, Ly),
                       z = (-Lz, 0),
                       topology=(Flat, Bounded, Bounded))

deltaN = 500kilometers   # North downwelling width
deltaS = 3000kilometers  # South upwelling width
deltaZ = 2000     # Vertical asymmetry
Ψᵢ(y, z)  = - 7 * ((1 - exp(-y / deltaS)) * (1 - exp(-(Ly - y) / deltaN)) * 
            sinpi(z / Lz) * exp(z/deltaZ))

Ψ = Field{Center, Face, Face}(grid)
set!(Ψ, Ψᵢ)
fill_halo_regions!(Ψ, arch)

# Set velocity field from streamfunction
v = YFaceField(grid)
w = ZFaceField(grid)
v .= - ∂z(Ψ)
w .= + ∂y(Ψ)
fill_halo_regions!(v, arch)
fill_halo_regions!(w, arch)

# (I have to specify u to allow CheckPointer)
u = XFaceField(grid) 
fill_halo_regions!(u, arch)

############################# Model setup ############################# 

kz(y,z,t) = 1e-4 + 5e-3 * (tanh((z+100)/20)+1) + 1e-2 * exp(-(z+4000)/50)
tracer_vertical_closure = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), κ=kz)
tracer_horizontal_closure = HorizontalScalarDiffusivity(κ=1e3)

# Model
model = HydrostaticFreeSurfaceModel(grid = grid,
                                    biogeochemistry = NutrientsPlanktonBacteriaDetritus(; grid),
                                    velocities = PrescribedVelocityFields(; u, v, w),
                                    tracers = (:N, :P, :Z, :B, :D1, :D2),
                                    coriolis = nothing,
                                    buoyancy = nothing,
                                    closure = (tracer_vertical_closure, tracer_horizontal_closure))

set!(model, N=3, P=1e-1, Z=1e-1, B=1e-1, D1=8e-2, D2=2e-2) 

# spinup_time = 365.25*2000days
# compute_time = 365days

simulation = Simulation(model; Δt = 3hour, stop_time=3days) 

# Print the progress 
progress(sim) = @printf("Iteration: %d, time: %s , total(N): %.2e\n",
            iteration(sim), prettytime(sim),
            sum(model.tracers.N) + sum(model.tracers.P) + sum(model.tracers.B) + sum(model.tracers.D1) + sum(model.tracers.D2))
add_callback!(simulation, progress, IterationInterval(100))

# outputs = (
#             N = model.tracers.N,
#             P = model.tracers.P,
#             Z = model.tracers.Z,
#             B = model.tracers.B,
#             D1 = model.tracers.D1,
#             D2 = model.tracers.D2
#             )

filename = "test1.jld2"
simulation.output_writers[:simple_output] =
        JLD2OutputWriter(model, model.tracers; 
                        filename,
                        schedule = TimeInterval(3days), 
                        overwrite_existing = true)

# simulation.output_writers[:checkpointer] = Checkpointer(model,
#             schedule = TimeInterval(365.25*500days),
#             prefix = "AMOC_NPZBD_checkpoint",
#             overwrite_existing = false)
        
run!(simulation) #, pickup = false)

###################################################################
########################## Visualization ##########################
###################################################################
#=
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
y = ynodes(Pt)

fig = Figure(; size = (1200, 800), fontsize = 18)

axN  = Axis(fig[1, 1], xlabel = "y (×10³ km)", ylabel = "z (m)", title="[Nutrient] (mmol m⁻³)",aspect=1)
axP  = Axis(fig[1, 3], xlabel="y (×10³ km)", ylabel = "z (m)", title="[Phytoplankton] (mmol m⁻³)",aspect=1)
axZ  = Axis(fig[1, 5], xlabel="y (×10³ km)", ylabel = "z (m)", title="[Zooplankton] (mmol m⁻³)",aspect=1)
axB  = Axis(fig[2, 1], xlabel="y (×10³ km)", ylabel = "z (m)", title="[Bacteria] (mmol m⁻³)",aspect=1)
axD1 = Axis(fig[2, 3], xlabel="y (×10³ km)", ylabel = "z (m)", title="[Dissolved Detritus] (mmol m⁻³)",aspect=1)
axD2 = Axis(fig[2, 5], xlabel="y (×10³ km)", ylabel = "z (m)", title="[Particulate Detritus] (mmol m⁻³)",aspect=1)

slider = Slider(fig[3, 1:5], range=1:nt, startvalue=1)
n = slider.value

title = @lift @sprintf("t = %d years", (t[$n] / 365day))
Label(fig[0, 1:5], title)

Nn  = @lift interior(Nt[$n], 1, :, :)
Pn  = @lift interior(Pt[$n], 1, :, :)
Zn  = @lift interior(Zt[$n], 1, :, :)
Bn  = @lift interior(Bt[$n], 1, :, :)
D1n = @lift interior(D1t[$n], 1, :, :)
D2n = @lift interior(D2t[$n], 1, :, :)

hm_N = heatmap!(axN, y./1e6, z, Nn; colorrange = (0,4),colormap = :thermal) 
Colorbar(fig[1, 2], hm_N; flipaxis = false)
ylims!(axN,-500,0)
hm_P = heatmap!(axP, y./1e6, z, Pn; colorrange = (0,0.2),colormap = :thermal) 
Colorbar(fig[1, 4], hm_P; flipaxis = false)
ylims!(axP,-500,0)
hm_Z = heatmap!(axZ, y./1e6, z, Zn; colorrange = (0,0.2),colormap = :thermal) 
# Colorbar(fig[1, 5], hm_Z; flipaxis = false)
ylims!(axZ,-500,0)
hm_B = heatmap!(axB, y./1e6, z, Bn; colorrange = (0,0.2),colormap = :thermal) 
Colorbar(fig[2, 2], hm_B; flipaxis = false)
ylims!(axB,-500,0)
hm_D1 = heatmap!(axD1, y./1e6, z, D1n; colorrange = (0,0.05),colormap = :thermal) 
Colorbar(fig[2, 4], hm_D1; flipaxis = false)
ylims!(axD1,-500,0)
hm_D2 = heatmap!(axD2, y./1e6, z, D2n; colorrange = (0,0.05),colormap = :thermal) 
# Colorbar(fig[2, 4], hm_D2; flipaxis = false)
ylims!(axD2,-500,0)

record(fig, "AMOC_NPZBD_50y.mp4", 1:nt, framerate=20) do nn
    n[] = nn
end
nothing #hide
=#
