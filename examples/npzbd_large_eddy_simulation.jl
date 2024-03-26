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

# Physical parameters
ρₒ = 1024 # kg m⁻³, reference density
cₚ = 3991 # heat capacity
T₀ = 20   # ᵒC, surface temperature
S₀ = 35   # g kg⁻¹, surface salinity
N² = 1e-5 # buoyancy frequency
g = 9.81 #Oceananigans.Buoyancy.g_Earth

# Boundary fluxes
heat_flux              = Q = 200 # W m⁻² (positive => cooling)
zonal_wind_stress      = τˣ = 0.0 # N m⁻²
meridional_wind_stress = τʸ = 0.0 # N m⁻²

# Output
simulation_name = "npzbd_large_eddy_simulation"
output_interval = 2minutes

grid = RectilinearGrid(size = (Nx, Ny, Nz),
                       halo = (4, 4, 4),
                       x = (0, Lx),
                       y = (0, Ly),
                       z = (-Lz, 0),
                       #topology = (Periodic, Periodic, Bounded))
                       topology = (Periodic, Bounded, Bounded))

shelf_width = 50 # meters
step(y, y₀, d) = (1 + tanh((y - y₀) / d)) / 2
shelf(x, y) = -Ly/2 + step(y, Ly/2, shelf_width)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(shelf))

T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Q / (ρₒ * cₚ)))
u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(τˣ / ρₒ))
v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(τʸ / ρₒ))

equation_of_state = TEOS10EquationOfState(reference_density=ρₒ)
buoyancy = SeawaterBuoyancy(; equation_of_state, gravitational_acceleration=g)
biogeochemistry = NutrientsPlanktonBacteriaDetritus(grid)

model = NonhydrostaticModel(; grid, buoyancy, biogeochemistry,
                            tracers = (:T, :S),
                            timestepper = :RungeKutta3,
                            advection = WENO(order=5),
                            boundary_conditions = (u=u_bcs, v=v_bcs, T=T_bcs))


α = SeawaterPolynomials.thermal_expansion(T₀, S₀, 0, equation_of_state)
β = SeawaterPolynomials.haline_contraction(T₀, S₀, 0, equation_of_state)
dTdz = N² / (α * g)
Tᵢ(x, y, z) = T₀ + dTdz * z

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
set!(model, u=ϵ, v=ϵ, w=ϵ, T=T₀, N=Nᵢ, P=Pᵢ, Z=Zᵢ, B=Bᵢ, D=Dᵢ)

simulation = Simulation(model, Δt=1.0, stop_time=1hour)

function progress(sim)
    msg = string("Iter: ", iteration(sim), ", time: ", prettytime(sim))
    msg *= string(", Δt: ", prettytime(sim.Δt))
    @info msg
    return nothing
end

add_callback!(simulation, progress, IterationInterval(10))
conjure_time_step_wizard!(simulation, cfl=0.5, IterationInterval(5))

xz_filename = simulation_name * "_xz.jld2"
outputs = merge(model.tracers, model.velocities)
simulation.output_writers[:xz] = JLD2OutputWriter(model, outputs,
                                                  schedule = TimeInterval(output_interval),
                                                  indices = (:, 1, :),
                                                  filename = xz_filename,
                                                  overwrite_existing = true)

run!(simulation)

xz_filename = simulation_name * "_xz.jld2"
Tt = FieldTimeSeries(xz_filename, "T")
Nt = FieldTimeSeries(xz_filename, "N")
Pt = FieldTimeSeries(xz_filename, "P")
Zt = FieldTimeSeries(xz_filename, "Z")
Bt = FieldTimeSeries(xz_filename, "B")
Dt = FieldTimeSeries(xz_filename, "D")
wt = FieldTimeSeries(xz_filename, "w")

# T = model.tracers.T
# N = model.tracers.N
# P = model.tracers.P
# Z = model.tracers.Z
# B = model.tracers.B
# D = model.tracers.D
# w = model.velocities.w

fig = Figure(size=(1200, 600))

axT = Axis(fig[1, 1], title="Temperature")
axw = Axis(fig[1, 2], title="Vertical velocity")
axN = Axis(fig[1, 3], title="Nutrients")
axP = Axis(fig[2, 1], title="Phytoplankton")
axZ = Axis(fig[2, 2], title="Zooplankton")
axB = Axis(fig[2, 3], title="Bacteria")
axD = Axis(fig[2, 4], title="Detritus")

Ntimes = length(Tt)

slider = Slider(fig[3, 1:4], range=1:Ntimes, startvalue=Ntimes)
n = slider.value

Tn = @lift interior(Tt[$n], :, 1, :)
wn = @lift interior(wt[$n], :, 1, :)
Nn = @lift interior(Nt[$n], :, 1, :)
Pn = @lift interior(Pt[$n], :, 1, :)
Zn = @lift interior(Zt[$n], :, 1, :)
Bn = @lift interior(Bt[$n], :, 1, :)
Dn = @lift interior(Dt[$n], :, 1, :)

heatmap!(axT, Tn)
heatmap!(axw, wn)
heatmap!(axN, Nn)
heatmap!(axP, Pn)
heatmap!(axZ, Zn)
heatmap!(axB, Bn)
heatmap!(axD, Dn)

display(fig)

