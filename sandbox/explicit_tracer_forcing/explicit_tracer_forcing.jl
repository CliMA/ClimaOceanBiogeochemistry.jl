using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using GLMakie

N² = 1e-5
Nz = 64
Lz = 256
f₀ = 1e-4
Jᵇ = 1e-8
Jᵘ = -1e-4
Nᴺ = 4
Nᴾ = 4

kᴺ = Tuple(1e-1 + 1e-1 * rand() for n = 1:Nᴺ)
μᵖ = Tuple(1e-1 + 1e-1 * rand() for n = 1:Nᴺ)

grid = RectilinearGrid(size=Nz, z=(-Lz, 0), topology=(Flat, Flat, Bounded))

top_b_bcs = FluxBoundaryCondition(Jᵇ)
top_u_bcs = FluxBoundaryCondition(Jᵘ)
b_bcs = FieldBoundaryConditions(top=top_b_bcs)
u_bcs = FieldBoundaryConditions(top=top_u_bcs)

nutrient_tracers = Tuple(Symbol(:N, n) for n = 1:Nᴺ)
plankton_tracers = Tuple(Symbol(:P, n) for n = 1:Nᴾ)

tracers = (:b, :e, nutrient_tracers..., plankton_tracers...)

nutrient_forcing_arrays = Tuple(zeros(size(grid)...) for n = 1:Nᴺ)
plankton_forcing_arrays = Tuple(zeros(size(grid)...) for n = 1:Nᴾ)

@inline array_forcing(i, j, k, grid, clock, fields, a) = @inbounds a[i, j, k]

nutrient_forcings = NamedTuple(name => Forcing(array_forcing, discrete_form=true, parameters=nutrient_forcing_arrays[n])
                               for (n, name) in enumerate(nutrient_tracers))
plankton_forcings = NamedTuple(name => Forcing(array_forcing, discrete_form=true, parameters=plankton_forcing_arrays[n])
                               for (n, name) in enumerate(plankton_tracers))

forcing = merge(nutrient_forcings, plankton_forcings)

model = HydrostaticFreeSurfaceModel(; grid, forcing, tracers,
                                    boundary_conditions = (; b=b_bcs, u=u_bcs),
                                    buoyancy = BuoyancyTracer(),
                                    coriolis = FPlane(f=f₀),
                                    closure = CATKEVerticalDiffusivity())

bᵢ(z) = N² * z

nutrient_initial_conditions = NamedTuple(name => 1e+1 for name in nutrient_tracers)
plankton_initial_conditions = NamedTuple(name => 1e-1 for name in plankton_tracers)

set!(model; b=bᵢ, nutrient_initial_conditions..., plankton_initial_conditions...)

simulation = Simulation(model, Δt=5minute, stop_time=1day)

function progress(sim)
    msg = string("Iter: ", iteration(sim), ", time: ", prettytime(sim))
    @info msg
end

add_callback!(simulation, progress, IterationInterval(10))

function update_forcing!(sim)
    # Do stuff
    return nothing
end

add_callback!(simulation, update_forcing!)

run!(simulation)

fig = Figure()
axu = Axis(fig[1, 1], xlabel="u (m s⁻¹)", ylabel="z (m)")
axb = Axis(fig[1, 2], xlabel="b (m s⁻²)", ylabel="z (m)")
axe = Axis(fig[1, 3], xlabel="e (m² s⁻²)", ylabel="z (m)")

u = model.velocities.u
b = model.tracers.b
e = model.tracers.e
z = znodes(u)

un = interior(u, 1, 1, :)
bn = interior(b, 1, 1, :)
en = interior(e, 1, 1, :)

lines!(axu, un, z)
lines!(axb, bn, z)
lines!(axe, en, z)



