using Oceananigans
using Oceananigans.Units
using ClimaOceanBiogeochemistry: NutrientsPlanktonBacteriaDetritus
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
simulation_name = "npzbd_large_eddy_simulation"
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

biogeochemistry = NutrientsPlanktonBacteriaDetritus(grid;)

model = NonhydrostaticModel(; grid, buoyancy, biogeochemistry,
                            tracers = (:T, :S),
                            timestepper = :RungeKutta3,
                            advection = WENO(order=5),
                            boundary_conditions = (u = u_bcs, T = T_bcs))


g = 9.81 #Oceananigans.Buoyancy.g_Earth
α = SeawaterPolynomials.thermal_expansion(T₀, S₀, 0, equation_of_state)

P₀ = 0.1 # μM
Z₀ = 0.1 # μM
Bᵢ = 0.1 # μM
Dᵢ = 0.1 # μM
N₀ = 0.001  # μM, surface nutrient concentration
hN = 10     # nutrient decay scale

Nᵢ(x, y, z) = N₀ * (1 - exp(z / hN))
Pᵢ(x, y, z) = P₀ *  exp(z / hN)
Zᵢ(x, y, z) = Z₀ *  exp(z / hN)

ϵ(x, y, z) = 1e-3 * randn()
set!(model, u=ϵ, v=ϵ, w=ϵ, T=T₀, N=Nᵢ, P=Pᵢ,Z=Zᵢ, B=Bᵢ, D=Dᵢ)

simulation = Simulation(model, Δt=1.0, stop_time=2hour) #, stop_iteration=100)

wizard = TimeStepWizard(cfl=0.5, max_change=1.1)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(5))

progress(sim) = @info string("Iter: ", iteration(sim), ", time: ", prettytime(sim))
simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

run!(simulation)

T = model.tracers.T
N = model.tracers.N
P = model.tracers.P
Z = model.tracers.Z
B = model.tracers.B
D = model.tracers.D
w = model.velocities.w

Tn = interior(T, :, 1, :)
wn = interior(w, :, 1, :)
Nn = interior(N, :, 1, :)
Pn = interior(P, :, 1, :)
Zn = interior(Z, :, 1, :)
Bn = interior(B, :, 1, :)
Dn = interior(D, :, 1, :)

fig = Figure()

axT = Axis(fig[1, 1])
axw = Axis(fig[1, 2])
axN = Axis(fig[1, 3])
axP = Axis(fig[2, 1])
axZ = Axis(fig[2, 2])
axB = Axis(fig[2, 3])
axD = Axis(fig[2, 4])

heatmap!(axT, Tn)
heatmap!(axw, wn)
heatmap!(axN, Nn)
heatmap!(axP, Pn)
heatmap!(axZ, Zn)
heatmap!(axB, Bn)
heatmap!(axD, Dn)

display(fig)

