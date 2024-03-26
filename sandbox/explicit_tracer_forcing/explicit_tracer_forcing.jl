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
Nᴺ = 2 # number of nutrient types
Nᴾ = 2 # number of phytoplankton types

kᴺ = Tuple(1e-1 + 1e-1 * rand() for n = 1:Nᴺ) # Range: 0.1 - 0.2
μᵖ = Tuple(1.0 + 1.0 * rand() for n = 1:Nᴺ) # Range: 1.0 - 2.0

grid = RectilinearGrid(size=Nz, z=(-Lz, 0), topology=(Flat, Flat, Bounded))

top_b_bcs = FluxBoundaryCondition(Jᵇ)
top_u_bcs = FluxBoundaryCondition(Jᵘ)
b_bcs = FieldBoundaryConditions(top=top_b_bcs)
u_bcs = FieldBoundaryConditions(top=top_u_bcs)

nutrient_tracers = Tuple(Symbol(:N, n) for n = 1:Nᴺ) # :N1, :N2, ...
plankton_tracers = Tuple(Symbol(:P, n) for n = 1:Nᴾ) # :P1, :P2, ...

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

#bᵢ(z) = N² * z

nutrient_initial_conditions = NamedTuple(name => 2e0 for name in nutrient_tracers)
plankton_initial_conditions = NamedTuple(name => 1e-1 for name in plankton_tracers)

#b=bᵢ,
set!(model; nutrient_initial_conditions..., plankton_initial_conditions...)

simulation = Simulation(model, Δt=5minute, stop_time=1day)

function progress(sim)
    msg = string("Iter: ", iteration(sim), ", time: ", prettytime(sim))
    @info msg
end

# TODO: add_callback: NOT DEFINED
# updated Oceananigans package on 03/25/2024 -> problem unresolved
add_callback!(simulation, progress, IterationInterval(10)) 

#=
∂t N1 = N1_rhs = -μᵖ[1]*(N1/(N1+kᴺ[1]))*P1 + mlin*P1
∂t N2 = N2_rhs

∂t P1 = P1_rhs
∂t P2 = P2_rhs

=#
mlin = 0.01

function update_forcing!(sim)
    Nx, Ny, Nz = size(sim.model.grid)
    N1_rhs = nutrient_forcing_arrays[1]
    N2_rhs = nutrient_forcing_arrays[2]
    P1_rhs = plankton_forcing_arrays[1]
    P2_rhs = plankton_forcing_arrays[2]
    for i = 1:Nx, j = 1:Ny, k = 1:Nz
        @inbounds begin
            N1_rhs[i, j, k] = -μᵖ[1]*(N1/(N1+kᴺ[1]))*P1 + mlin*P1
            N2_rhs[i, j, k] = -μᵖ[2]*(N2/(N2+kᴺ[2]))*P2 + mlin*P2
            P1_rhs[i, j, k] = μᵖ[1]*(N1/(N1+kᴺ[1]))*P1 - mlin*P1
            P2_rhs[i, j, k] = μᵖ[2]*(N2/(N2+kᴺ[2]))*P2 - mlin*P2
        end
    end
    # Do stuff
    N1 += N1_rhs
    return nothing
end

add_callback!(simulation, update_forcing!)

run!(simulation)

fig = Figure()
#axu = Axis(fig[1, 1], xlabel="u (m s⁻¹)", ylabel="z (m)")
#axb = Axis(fig[1, 2], xlabel="b (m s⁻²)", ylabel="z (m)")
#axe = Axis(fig[1, 3], xlabel="e (m² s⁻²)", ylabel="z (m)")
axN = Axis(fig[1, 1], xlabel="N (mmol m⁻³)", ylabel="z (m)")
axP = Axis(fig[1, 2], xlabel="P (mmol m⁻³)", ylabel="z (m)")

u = model.velocities.u # horizontal velocity (flow rate)
b = model.tracers.b # buoyancy
e = model.tracers.e # turbulence kinetic energy

N1 = model.tracers.N1
N2 = model.tracers.N2
P1 = model.tracers.P1
P2 = model.tracers.P2
z = znodes(N1)

#un = interior(u, 1, 1, :)
#bn = interior(b, 1, 1, :)
#en = interior(e, 1, 1, :)
N1n = interior(N1, 1, 1, :)
N2n = interior(N2, 1, 1, :)
P1n = interior(P1, 1, 1, :)
P2n = interior(P2, 1, 1, :)

lines!(axN, N1n, z)
lines!(axN, N2n, z)

lines!(axP, P1n, z)
lines!(axP, P2n, z)

save("test.png", fig)



