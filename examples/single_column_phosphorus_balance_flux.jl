# The phosphorus cycle: PO₄, DOP, POP
#
# This example illustrates how to use ClimaOceanBiogeochemistry's
# `CarbonAlkalinityNutrients` model in a single column context.

using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients

using Oceananigans
using Oceananigans.Units
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Advection: _advective_tracer_flux_z
using Oceananigans.Grids: znode
using Oceananigans.Forcings: maybe_constant_field
using Printf
using CairoMakie

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

@inline κ(z, t) = 1e-4 + 1e-2 * exp(z / 50) + 1e-2 * exp(-(z + 1000) / 50)
vertical_diffusion = VerticalScalarDiffusivity(; κ)

# The following two lines should be later debugged in Oceananigans
import Oceananigans.Models.HydrostaticFreeSurfaceModels: validate_tracer_advection, AbstractAdvectionScheme, SingleColumnGrid
validate_tracer_advection(tracer_advection::AbstractAdvectionScheme, ::SingleColumnGrid) = tracer_advection, NamedTuple()

# We put the pieces together.
# The important line here is `biogeochemistry = CarbonAlkalinityNutrients()`.
# We use all default parameters.
 
wᵢ = -10/day
w = ZFaceField(grid)
set!(w, wᵢ)
# Manually set boundary conditions
w[1, 1, 1] = 0
w[1, 1, Nz+1] = 0

# Function: the flux of the sinking tracer b
@inline function a_flux(i, j, grid, clock, fields, parameters)
    advection = parameters.advection
    k_bound = parameters.k_bound
    w_sinking = parameters.w_sinking
    wb = _advective_tracer_flux_z(i, j, k_bound, grid, advection, w_sinking, fields.POP)
    return wb*0.1
end

# k_bound: the depth at which POP flux sinks out (80 = 800 m)
parameters = (; advection = WENO(), w_sinking=w, k_bound=100)

# Add the b flux out as a top boundary flux of a
a_flux_bc = FluxBoundaryCondition(a_flux, discrete_form=true; parameters)
top_flux_bcs = FieldBoundaryConditions(top=a_flux_bc)
bottom_flux_bcs = FieldBoundaryConditions(bottom=a_flux_bc)

T_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(8e-4))
#                                bottom = GradientBoundaryCondition(0.0001))

# A pulsed flux as a function of time
# Fᵖ(t) = ifelse(t % (30days) <(3days), -1e-6, 0.0)
# top_pulse_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Fᵖ))

model = HydrostaticFreeSurfaceModel(; grid,
                                    velocities = PrescribedVelocityFields(),
                                    biogeochemistry = CarbonAlkalinityNutrients(; grid,
                                                                                #maximum_net_community_production_rate  = 5e-2/day,
                                                                                fraction_of_particulate_export  = 1,
                                                                                #dissolved_organic_phosphorus_remin_timescale  = 3/365.25days,
                                                                                #martin_curve_exponent = 0.84,
                                                                                particulate_organic_phosphorus_sinking_velocity = -10/day),
                                    tracers = (:DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe),
                                    #boundary_conditions = (; POP=bottom_flux_bcs, PO₄=top_flux_bcs),#POP=bottom_flux_bcs, PO₄=top_flux_bcs), 
                                    tracer_advection = WENO(),
                                    buoyancy = nothing,
                                    closure = vertical_diffusion) 

# ## Initial conditions
#
# We initialize the model with reasonable concentrations

N₀ = 1e-2 # Surface nutrient concentration
D₀ = 1e-3 # Surface detritus concentration

Nᵢ(z) = N₀ * (-z / 115).^0.84
Dᵢ(z) = D₀ * (-z / 115).^-0.84

# set!(model, DIC=2.1, ALK=2.35, NO₃=3.6e-2, PO₄=2e-3, DOP=5e-4, POP=1.5e-3, Fe = 1e-6) # mol PO₄ m⁻³
set!(model, DIC=2.1, ALK=2.35, NO₃=1e-2, PO₄=Nᵢ, DOP=0, POP=Dᵢ, Fe = 1e-5) # mol PO₄ m⁻³

simulation = Simulation(model, Δt=10minutes, stop_time=365.25days)

function progress(sim)
    @printf("Iteration: %d, time: %s, total(P): %.2e \n",
            iteration(sim), prettytime(sim),
            sum(model.tracers.PO₄)+sum(model.tracers.DOP)+sum(model.tracers.POP))
    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

# Add a callback to calculate rates (production or remineralization) at each timestep

filename = "po4_dop_pop.jld2"

simulation.output_writers[:fields] = JLD2OutputWriter(model, model.tracers;
                                                      filename,
                                                      schedule = TimeInterval(1day),
                                                      overwrite_existing = true)

run!(simulation)

# ###################### Visualization ######################

# All that's left is to visualize the results.
PO4t = FieldTimeSeries(filename, "PO₄")
DOPt = FieldTimeSeries(filename, "DOP")
POPt = FieldTimeSeries(filename, "POP")

t = PO4t.times
nt = length(t)
z = znodes(PO4t)

fig = Figure(;size=(800, 600))#resolution=(1200, 600))

ax1 = Axis(fig[1, 1], ylabel="z (m)", xlabel="[Inorganic PO₄] (mol m⁻³)")
ax2 = Axis(fig[1, 2], ylabel="z (m)", xlabel="[Organic P] (mol m⁻³)")
# ax3 = Axis(fig[1, 3], ylabel="z (m)", xlabel="ln[POP] at end point")

xlims!(ax1, -0.0001, 0.06)
xlims!(ax2, 0, 0.001)
# xlims!(ax3, -20, -7)
ylims!(ax1, -1000, 0)
ylims!(ax2, -1000, 0)
# ylims!(ax3, -1000, 0)

slider = Slider(fig[2, 1:2], range=1:nt, startvalue=1)
n = slider.value

title = @lift @sprintf("Equilibrium biogeochemistry at t = %d days", t[$n] / day)
Label(fig[0, 1:2], title)

PO4n = @lift interior(PO4t[$n], 1, 1, :)
DOPn = @lift interior(DOPt[$n], 1, 1, :)
POPn = @lift interior(POPt[$n], 1, 1, :)

# POP_ln = log.(interior(POPt[end], 1, 1, :))

lines!(ax1, PO4n, z, label = "PO₄")
axislegend(ax1, position = :rt)
# lines!(ax2, DOPn, z, label = "DOP")
lines!(ax2, POPn, z, label = "POP")
axislegend(ax2, position = :rb)
# lines!(ax3, POP_ln, z, label = "ln(POP)")
# # lines!(ax3, 0.84/z.-10, z, label = "-b/z")
# axislegend(ax3, position = :lb)

record(fig, "po4_dop_pop.mp4", 1:nt, framerate=24) do nn
    n[] = nn
end
nothing

#=
POP_last = interior(POPt[end], 1, 1, :) 
POC_last = POP_last * 106
POC_flux = POC_last * 5 # flux = concentration x sinking rate

last_frame = Figure(size=(600, 600))

axF  = Axis(last_frame[1, 1], xlabel="POC flux (mmol m⁻² d⁻¹)", ylabel="Depth (m)")
xlims!(axF, 0, 0.01)
ylims!(axF, -800, -100)
lines!(axF, POC_flux, z, linewidth = 3,label = "Modeled flux")

z₀ = z[86]# z[91] = 95 m; z[90] = 105 m
martin = POC_flux[86].*(z./z₀).^-0.84

lines!(axF, martin, z, linewidth = 3,label = "Martin curve")
axislegend(position = :rb)
display(last_frame)
save("po4_dop_pop.png", last_frame)
nothing #hide
=#