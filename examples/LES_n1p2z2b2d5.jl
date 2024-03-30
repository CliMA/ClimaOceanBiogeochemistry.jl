using Oceananigans
using Oceananigans.Units
using ClimaOceanBiogeochemistry: MultiNPZBD
using SeawaterPolynomials
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using GLMakie

# Resolution
Nx = Ny = 64
Nz = 64

# Domain
Lx = 100
Ly = 200 # meters
Lz = 100 # meters

# Output
simulation_name = "LES_multinpzbd"
output_interval = 2minutes

grid = RectilinearGrid(size = (Nx, Ny, Nz),
                       x = (0, Lx),
                       y = (0, Ly),
                       z = (-Lz, 0),
                       topology = (Periodic, Periodic, Bounded))


ρₒ = 1024
cₚ = 3991
T₀ = 20 # ᵒC, surface temperature
S₀ = 35 # psu, surface salinity

equation_of_state = TEOS10EquationOfState(reference_density=ρₒ)

# Boundary fluxes
heat_flux = Q = 200 # W m⁻²
wind_stress = τˣ = 0.0 # N m⁻²

T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Q / (ρₒ * cₚ)))
u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(τˣ / ρₒ))

buoyancy = SeawaterBuoyancy(; equation_of_state)

biogeochemistry = MultiNPZBD(grid;)

model = NonhydrostaticModel(; grid, buoyancy, biogeochemistry,
                            tracers = (:T, :S,:N1,:P1,:P2,:Z1,:Z2,:B1,:B2,:D1,:D2,:D3,:D4,:D5),
                            timestepper = :RungeKutta3,
                            advection = WENO(order=5),
                            boundary_conditions = (u = u_bcs, T = T_bcs))


g = 9.81 #Oceananigans.Buoyancy.g_Earth
α = SeawaterPolynomials.thermal_expansion(T₀, S₀, 0, equation_of_state)

P₀ = 1.0 # μM
Z₀ = 0.8 # μM
Bᵢ = 0.1 # μM
D₀ = 0.5 # μM
N₀ = 0.001  # μM, surface nutrient concentration
hN = 10     # nutrient decay scale

Nᵢ(x, y, z) = N₀ * (1 - exp(z / hN))
Pᵢ(x, y, z) = P₀ * (1 - exp(z/10))
Zᵢ(x, y, z) = Z₀ * (1 - exp(z/10))
Dᵢ(x, y, z) = D₀ * (1 - exp(z/10))

ϵ(x, y, z) = 1e-3 * randn()

set!(model, u=ϵ, v=ϵ, w=ϵ, T=20, 
    N1=Nᵢ, P1=Pᵢ,P2=Pᵢ,
    Z1=Zᵢ,Z2=Zᵢ,
    B1=1.5.*Bᵢ,B2=Bᵢ, 
    D1=Dᵢ, D2=Dᵢ, D3=Dᵢ, D4=Dᵢ, D5=Dᵢ)

simulation = Simulation(model, Δt=1.0, stop_time=1hour) #stop_iteration=100)

wizard = TimeStepWizard(cfl=0.5, max_change=1.1)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(5))

progress(sim) = @info string("Iter: ", iteration(sim), ", time: ", prettytime(sim))
simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

xz_filename = simulation_name * "_xz.jld2"
outputs = merge(model.tracers, model.velocities)
simulation.output_writers[:xz] = JLD2OutputWriter(model, outputs,
                                                  schedule = TimeInterval(output_interval),
                                                  indices = (:, 1, :),
                                                  filename = xz_filename,
                                                  overwrite_existing = true)

run!(simulation)

T = model.tracers.T
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

w = model.velocities.w

Tn = interior(T, :, 1, :)
wn = interior(w, :, 1, :)
N1n = interior(N1, :, 1, :)
P1n = interior(P1, :, 1, :)
P2n = interior(P2, :, 1, :)
Z1n = interior(Z1, :, 1, :)
Z2n = interior(Z2, :, 1, :)
B1n = interior(B1, :, 1, :)
B2n = interior(B2, :, 1, :)
D1n = interior(D1, :, 1, :)
D2n = interior(D2, :, 1, :)
D3n = interior(D3, :, 1, :)
D4n = interior(D4, :, 1, :)
D5n = interior(D5, :, 1, :)

fig = Figure()

axT = Axis(fig[1, 1], title="Temperature")
axw = Axis(fig[1, 2], title="Vertical velocity")
axN1 = Axis(fig[1, 3], title="Nutrients")
axP1 = Axis(fig[1, 4], title="P_slow")
axP2 = Axis(fig[1, 5], title="P_fast")
axZ1 = Axis(fig[1, 6], title="Z_P")
axZ2 = Axis(fig[1, 7], title="Z_B")
axB1 = Axis(fig[2, 1], title="B_slow")
axB2 = Axis(fig[2, 2], title="B_fast")
axD1 = Axis(fig[2, 3], title="LDOM")
axD2 = Axis(fig[2, 4], title="SLDOM")
axD3 = Axis(fig[2, 5], title="RDOM")
axD4 = Axis(fig[2, 6], title="POM_slow")
axD5 = Axis(fig[2, 7], title="POM_fast")

heatmap!(axT, Tn)
heatmap!(axw, wn)
heatmap!(axN1, N1n)
heatmap!(axP1, P1n)
heatmap!(axP2, P2n)
heatmap!(axZ1, Z1n)
heatmap!(axZ2, Z2n)
heatmap!(axB1, B1n)
heatmap!(axB2, B2n)
heatmap!(axD1, D1n)
heatmap!(axD2, D2n)
heatmap!(axD3, D3n)
heatmap!(axD4, D4n)
heatmap!(axD5, D5n)

display(fig)

