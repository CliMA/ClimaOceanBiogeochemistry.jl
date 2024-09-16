# # test a 2D physical setting with Lagrangian particles

using GLMakie
using Printf
using Statistics

using Oceananigans
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Oceananigans.Units

#################################### Grid ####################################
# Parameters
Nx = 100 
Nz = 200
Lx = 1000kilometers   # m
Lz = 2000           # m

# We use a two-dimensional grid, with a `Flat` `y`-direction:
grid = RectilinearGrid(size = (Nx, Nz),
                       x = (0, Lx),
                       z = (-Lz, 0),
                       topology=(Bounded, Flat, Bounded))

#################################### Boundary conditions ####################################

# Wind stress boundary condition
# τ₀ = 1e-5           # m² s⁻²
@inline function τy(x, t, p)
    return 0.0125/1e3 * (sinpi(x/p.Lx/0.8-0.5)+1)*(-tanh(x/p.Lx*pi/0.1-9*pi)+1)
end

y_wind_stress = FluxBoundaryCondition(τy, parameters=(; Lx)) #τ₀,
v_bcs = FieldBoundaryConditions(top=y_wind_stress)

#################################### Model ####################################

background_diffusivity = ScalarDiffusivity(ν=1e-4, κ=1e-4)
catke = CATKEVerticalDiffusivity()
# horizontal_closure = HorizontalScalarDiffusivity(ν=1e3, κ=1e3)

# Add Lagrangian Particles
n_particles = 3
x₀ = Lx/2 * ones(n_particles) #* collect([0.1,0.5,0.9]) #rand(n_particles)
y₀ = ones(n_particles)
z₀ = -Lz * collect([0.005, 0.05, 0.1]) #* ones(n_particles)

# using StructArrays
# struct LagrangianP{T}
#     x :: T
#     y :: T
#     z :: T
#     # Particles_2D :: P
# end
# Particles_init = ones(n_particles)
# particles = StructArray{LagrangianP}((x₀, y₀, z₀))
# lagrangian_particles = LagrangianParticles(particles; tracked_fields=(; Particles_2D))
lagrangian_particles = LagrangianParticles(x=x₀, y=y₀, z=z₀) 

# Model
model = HydrostaticFreeSurfaceModel(; grid,
                                    buoyancy = BuoyancyTracer(),
                                    coriolis = FPlane(; f=1e-4),
                                    closure = (catke, background_diffusivity), #horizontal_closure),
                                    tracers = (:b, :e),
                                    momentum_advection = WENO(),
                                    particles = lagrangian_particles,
                                    boundary_conditions = (; v=v_bcs))

M² = 0         # 1e-7 s⁻², squared buoyancy frequency
N² = 0         # 1e-5 s⁻²
bᵢ(x, z) = N² * z + M² * x

set!(model, b=bᵢ) 

simulation = Simulation(model; Δt = 30minutes, stop_time=365.25days)

# outputs = (u = model.velocities.u,
#             w = model.velocities.w,
#             c = model.tracers.c,
#             h = model.tracers.h,
#             p = model.particles)

# simulation.output_writers[:simple_output] =
#     JLD2OutputWriter(model, outputs, 
#                      schedule = TimeInterval(1days), #30minutes
#                      filename = "PassiveTracer_2D",
#                      overwrite_existing = true)

run!(simulation)

############################ Visualizing the solution ############################

xp = model.particles.properties.x
zp = model.particles.properties.z

fig = Figure(size = (500, 500))

ax_p = Axis(fig[1, 1]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
x₀ = Lx/2 * ones(n_particles) 
z₀ = collect([-10, -100, -200]) 
# scatter!(ax_p, x₀,z₀, color =:blue, markersize = 15, label="initial") 
# scatter!(ax_p, xp,zp, color =:red, markersize = 10, label="final") 
scatter!(ax_p, x₀[1],z₀[1], color =:blue, markersize = 15, label="P1 initial") 
scatter!(ax_p, xp[1],zp[1], color =:blue, markersize = 15, marker =:xcross,label="P1 final") 
scatter!(ax_p, x₀[2],z₀[2], color =:black, markersize = 15, label="P2 initial") 
scatter!(ax_p, xp[2],zp[2], color =:black, markersize = 15, marker =:xcross,label="P2 final") 
scatter!(ax_p, x₀[3],z₀[3], color =:red, markersize = 15, label="P3 initial") 
scatter!(ax_p, xp[3],zp[3], color =:red, markersize = 15, marker =:xcross,label="P3 final") 

# xlims!(ax_p, 0, Lx)
axislegend(ax_p, position = :lb)
display(fig)


