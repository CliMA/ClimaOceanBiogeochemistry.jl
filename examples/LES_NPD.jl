using Oceananigans
using Oceananigans.Units
using ClimaOceanBiogeochemistry: NutrientsPlanktonBacteriaDetritus
using SeawaterPolynomials
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using Statistics
using GLMakie

# Resolution
Nx = Ny = 32
Nz = 32

# Domain
Lx = 1000
Ly = 1000 # meters
Lz = 500 # meters

# Output
simulation_name = "LES_test"
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

biogeochemistry = NutrientsPlanktonBacteriaDetritus(; grid, 
                                                    linear_remineralization_rate = 0.15/day, 
                                                    # maximum_bacteria_growth_rate = 1.2/day,
                                                    # detritus_half_saturation     = 0.5,
                                                    # bacteria_yield               = 0.5,
                                                    detritus_vertical_velocity = -5/day)

model = NonhydrostaticModel(; grid, buoyancy, biogeochemistry,
                            tracers = (:T, :S),
                            timestepper = :RungeKutta3,
                            advection = WENO(order=5),
                            boundary_conditions = (u = u_bcs, T = T_bcs))


g = 9.81 #Oceananigans.Buoyancy.g_Earth
α = SeawaterPolynomials.thermal_expansion(T₀, S₀, 0, equation_of_state)

P₀ = 0.1 # μM
#Z₀ = 0 # μM
#Bᵢ = 0.1 # μM
D₀ = 0.5 # μM
N₀ = 1  # μM, surface nutrient concentration
hN = 10     # nutrient decay scale

Nᵢ(x, y, z) = N₀ * (1 - exp(z / hN))
Pᵢ(x, y, z) = P₀ *  exp(z / hN)
Dᵢ(x, y, z) = D₀ *  exp(z / hN)
#Zᵢ(x, y, z) = Z₀ *  exp(z / hN)

ϵ(x, y, z) = 1e-3 * randn()
set!(model, u=ϵ, v=ϵ, w=ϵ, T=T₀, N=Nᵢ, P=Pᵢ,Z=0, B=0, D=Dᵢ)

simulation = Simulation(model, Δt=1.0, stop_time=3hour) #, stop_iteration=100)

wizard = TimeStepWizard(cfl=0.5, max_change=1.1)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(5))

progress(sim) = @info string("Iter: ", iteration(sim), ", time: ", prettytime(sim))
simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

xz_filename = simulation_name * "_xz.jld2"
outputs = merge(model.tracers, model.velocities)
#=
saved_data = Dict(
    "model" => model,  # Include the model object in the saved data
    "outputs" => merge(model.tracers, model.velocities)
    )
    =#
simulation.output_writers[:xz] = JLD2OutputWriter(model, outputs,
                                                  schedule = TimeInterval(output_interval),
                                                  indices = (:, 1, :),
                                                  filename = xz_filename,
                                                  overwrite_existing = true)

run!(simulation)

#T = model.tracers.T
N = model.tracers.N
P = model.tracers.P
#Z = model.tracers.Z
#B = model.tracers.B
D = model.tracers.D
w = model.velocities.w

#Tn = interior(T, :, 1, :)
wn = interior(w, :, 1, :)
Nn = interior(N, :, 1, :)
Pn = interior(P, :, 1, :)
#Zn = interior(Z, :, 1, :)
#Bn = interior(B, :, 1, :)
Dn = interior(D, :, 1, :)

fig1 = Figure()

#axT = Axis(fig1[1, 1], title="Temperature")
axw = Axis(fig1[1, 1], title="Vertical velocity")
axN = Axis(fig1[1, 2], title="Nutrients")
axP = Axis(fig1[1, 3], title="Phytoplankton")
#axZ = Axis(fig[2, 2], title="Zooplankton")
#axB = Axis(fig1[1, 3], title="Bacteria")
axD = Axis(fig1[1, 4], title="Detritus")

#heatmap!(axT, Tn)
heatmap!(axw, wn)
heatmap!(axN, Nn)
heatmap!(axP, Pn)
#heatmap!(axZ, Zn)
#heatmap!(axB, Bn)
heatmap!(axD, Dn)

# Plot average profiles
w_col = vec(mean(wn, dims=1))
N_col = vec(mean(Nn, dims=1))
P_col = vec(mean(Pn, dims=1))
#B_col = vec(mean(Bn, dims=1))
D_col = vec(mean(Dn, dims=1))

axwcol = Axis(fig1[2, 1], title="Vertical velocity")
axNcol = Axis(fig1[2, 2], title="Nutrients")
axPcol = Axis(fig1[2, 3], title="Phytoplankton")
#axBcol = Axis(fig1[2, 3], title="Bacteria")
axDcol = Axis(fig1[2, 4], title="Detritus")

z = vec(reshape(collect(1.0:32.0), 1, 32))
lines!(axwcol, w_col[2:end], z)
lines!(axNcol, N_col, z)
lines!(axPcol, P_col, z)
#lines!(axBcol, B_col, z)
lines!(axDcol, D_col, z)

#display(fig1)
save("NPZDB_13.png", fig1)
