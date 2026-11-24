# # One-degree global ocean simulation
#
# This example configures a global ocean--sea ice simulation at 1ᵒ horizontal resolution with
# realistic bathymetry and some closures.
#
# For this example, we need Oceananigans, ClimaOcean, OrthogonalSphericalShellGrids, and
# CairoMakie to visualize the simulation. Also we need CFTime and Dates for date handling.

using ClimaOcean
using ClimaOcean.ECCO
using Oceananigans
using Oceananigans.Units
using OrthogonalSphericalShellGrids
using CFTime
using Dates
using Printf

using ClimaOceanBiogeochemistry
using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients
using ClimaOceanBiogeochemistry.CarbonSystemSolvers.UniversalRobustCarbonSolver: UniversalRobustCarbonSystem
using ClimaOceanBiogeochemistry.CarbonSystemSolvers: CarbonSystemParameters, CarbonSolverParameters, CarbonCoefficientParameters
using Oceananigans.Grids: znode

arch = GPU()

# ### Grid and Bathymetry

Nx = 360
Ny = 180
Nz = 100

r_faces = exponential_z_faces(; Nz, depth=5000, h=34)
z_faces = Oceananigans.MutableVerticalDiscretization(r_faces)

underlying_grid = TripolarGrid(arch;
                               size = (Nx, Ny, Nz),
                               z = z_faces,
                               halo = (5, 5, 4),
                               first_pole_longitude = 70,
                               north_poles_latitude = 55)

bottom_height = regrid_bathymetry(underlying_grid;
                                  minimum_depth = 10,
                                  interpolation_passes = 75,
                                  major_basins = 2)

# For this bathymetry at this horizontal resolution we need to manually open the Gibraltar strait.
view(bottom_height, 102:103, 124, 1) .= -400
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height); active_cells_map=true)

# ### Restoring
#
# We include temperature and salinity surface restoring to ECCO data.

restoring_rate  = 1 / 10days
z_below_surface = r_faces[end-1]

mask = LinearlyTaperedPolarMask(southern=(-80, -70), northern=(70, 90), z=(z_below_surface, 0))

dates = DateTimeProlepticGregorian(1993, 1, 1) : Month(1) : DateTimeProlepticGregorian(1994, 1, 1)
temperature = ECCOMetadata(:temperature; dates, version=ECCO4Monthly(), dir="./")
salinity    = ECCOMetadata(:salinity;    dates, version=ECCO4Monthly(), dir="./")

FT = ECCORestoring(temperature, grid; mask, rate=restoring_rate)
FS = ECCORestoring(salinity,    grid; mask, rate=restoring_rate)
#forcing = (T=FT, S=FS)

# ### Closures
# We include a Gent-McWilliam isopycnal diffusivity as a parameterization for the mesoscale
# eddy fluxes. For vertical mixing at the upper-ocean boundary layer we include the CATKE
# parameterization. We also include some explicit horizontal diffusivity.

using Oceananigans.TurbulenceClosures: IsopycnalSkewSymmetricDiffusivity,
                                       DiffusiveFormulation

eddy_closure = IsopycnalSkewSymmetricDiffusivity(κ_skew=1e3, κ_symmetric=1e3,
                                                 skew_flux_formulation=DiffusiveFormulation())
vertical_mixing = ClimaOcean.OceanSimulations.default_ocean_closure()

closure = (eddy_closure, vertical_mixing)

# biogeochemistry
#CO₂_flux = Field{Center, Center, Nothing}(grid)
iron = ECCOMetadata(:Fe;  dates, version=ECCO4DarwinMonthly(), dir="./")
imask= LinearlyTaperedPolarMask(southern=(-80, -80), northern=(90, 90), z=(z_below_surface, 0)) 
irate= 1/5days
FF   = ECCORestoring(iron, grid; mask=imask, rate=irate)

# Initialatize the carbon and alkalinity FW forcing (related to salinity restoring)
#FF = FieldBoundaryConditions(
#    top = FluxBoundaryCondition(DustDeposition)
#    )
FC  = Field{Center, Center, Nothing}(grid)
FA  = Field{Center, Center, Nothing}(grid)
FF2 = Field{Center, Center, Nothing}(grid)
forcing = (T=FT, S=FS) #, Fe=FF)

# Initialize carbon surface boundary condition (CO2 fluxes)
BC  = FieldBoundaryConditions(
    top = FluxBoundaryCondition(FC),
    )
BA = FieldBoundaryConditions(
    top = FluxBoundaryCondition(FA),
    )
#BF = FieldBoundaryConditions(
#    top = FluxBoundaryCondition(FF2),
#    )
# Need to specify iron boundary condition otherwise iron forcing has nowhere to go
boundary_conditions = (DIC=BC, ALK=BA) #, Fe=BF)

biogeochemistry = CarbonAlkalinityNutrients(; 
            grid,
            maximum_net_community_production_rate = 0.5/365.25days,
            carbon_solver_params = (Sᵒᵖᵗˢ = CarbonSolverParameters(),),
            atmospheric_pCO₂ = 380e-6,
)

# ### Ocean simulation
# Now we bring everything together to construct the ocean simulation.
# We use a split-explicit timestepping with 30 substeps for the barotropic
# mode.

free_surface = SplitExplicitFreeSurface(grid; substeps=30)

momentum_advection = WENOVectorInvariant(vorticity_order=3)
tracer_advection   = Centered()

ocean = ocean_simulation(grid;
                         momentum_advection,
                         tracer_advection,
                         closure,
                         forcing,
                         free_surface,
			 boundary_conditions,
			 biogeochemistry,
			 tracers = (:T, :S, :e, :DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe))

# ### Initial condition

# We initialize the ocean from the ECCO state estimate.

set!(ocean.model, T=ECCOMetadata(:temperature; dates=first(dates)),
                  S=ECCOMetadata(:salinity;    dates=first(dates)),
                DIC=ECCOMetadata(:DIC; dates=first(dates), version=ECCO4DarwinMonthly()),
                ALK=ECCOMetadata(:ALK; dates=first(dates), version=ECCO4DarwinMonthly()),
                PO₄=ECCOMetadata(:PO₄; dates=first(dates), version=ECCO4DarwinMonthly()),
                NO₃=ECCOMetadata(:NO₃; dates=first(dates), version=ECCO4DarwinMonthly()),
                DOP=ECCOMetadata(:DOP; dates=first(dates), version=ECCO4DarwinMonthly()),
                POP=ECCOMetadata(:POP; dates=first(dates), version=ECCO4DarwinMonthly()),
                Fe =ECCOMetadata(:Fe;  dates=first(dates), version=ECCO4DarwinMonthly()),
              )
# ### Atmospheric forcing

# We force the simulation with an JRA55-do atmospheric reanalysis.
radiation  = Radiation(arch)
atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(20))

# ### Coupled simulation

# Now we are ready to build the coupled ocean--sea ice model and bring everything
# together into a `simulation`.

# We use a relatively short time step initially and only run for a few days to
# avoid numerical instabilities from the initial "shock" of the adjustment of the
# flow fields.

coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)
simulation = Simulation(coupled_model; Δt=1minutes, stop_time=10days)

# Copy incident shortwave ratiation to PAR field and scaled salinity flux to DIC/ALK forcing
ClimaOceanBiogeochemistry.copy_surface_atmospheric_state_for_bgc!(simulation)

# Initialize pCO2 and PAR properly
ClimaOceanBiogeochemistry.update_biogeochemical_state!(
	simulation.model.ocean.model.biogeochemistry, 
	simulation.model.ocean.model,
)

# ### A progress messenger
#
# We write a function that prints out a helpful progress message while the simulation runs.

wall_time = Ref(time_ns())

function progress(sim)
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T
    Tmax = maximum(interior(T))
    Tmin = minimum(interior(T))
    umax = (maximum(abs, interior(u)),
            maximum(abs, interior(v)),
            maximum(abs, interior(w)))

    step_time = 1e-9 * (time_ns() - wall_time[])

    msg1 = @sprintf("time: %s, iteration: %d, Δt: %s, ", prettytime(sim), iteration(sim), prettytime(sim.Δt))
    msg2 = @sprintf("max|u|: (%.2e, %.2e, %.2e) m s⁻¹, ", umax...)
    msg3 = @sprintf("extrema(T): (%.2f, %.2f) ᵒC, ", Tmax, Tmin)
    msg4 = @sprintf("wall time: %s \n", prettytime(step_time))

    @info msg1 * msg2 * msg3 * msg4

     wall_time[] = time_ns()

     return nothing
end

# And add it as a callback to the simulation.
add_callback!(simulation, progress, IterationInterval(10))

# UpdateStateCallsite callbacks occur before BGC update state, where PAR and CO2 fluxes are updated.
#    THEN tendencies are computed (to use the PAR).
add_callback!(simulation, ClimaOceanBiogeochemistry.copy_surface_atmospheric_state_for_bgc!,
			  callsite   = Oceananigans.UpdateStateCallsite(), 
)

#Air-sea CO2 flux parameterization
include("calculate_air_sea_carbon_exchange.jl")
add_callback!(simulation, calculate_air_sea_carbon_exchange!, 
			  callsite   = Oceananigans.UpdateStateCallsite(), 
			 )

# ### Output
#
# We are almost there! We need to save some output. Below we choose to save _only surface_
# fields using the `indices` keyword argument. We save all velocity and tracer components.
# Note, that besides temperature and salinity, the CATKE vertical mixing parameterization
# also uses a prognostic turbulent kinetic energy, `e`, to diagnose the vertical mixing length.

outputs = merge(ocean.model.tracers, ocean.model.velocities,
       	     (CO₂_flux   = simulation.model.ocean.model.biogeochemistry.CO₂_flux, 
              ocean_pCO₂ = simulation.model.ocean.model.biogeochemistry.ocean_pCO₂, 
	      pH         = simulation.model.ocean.model.biogeochemistry.pH,
	      PAR        = simulation.model.ocean.model.biogeochemistry.PAR,
		),
)
ocean.output_writers[:surface] = JLD2OutputWriter(ocean.model, outputs;
                                                  schedule = TimeInterval(5days),
                                                  filename = "global_surface_fields_with_bgc",
                                                  indices = (:, :, grid.Nz),
                                                  with_halos = true,
                                                  overwrite_existing = true,
                                                  array_type = Array{Float32})

# ### Ready to run

# We are ready to press the big red button and run the simulation.

# After we run for a short time (here we set up the simulation with `stop_time = 10days`),
# we increase the timestep and run for longer.

run!(simulation)

simulation.Δt = 15minutes
simulation.stop_time = 360days

run!(simulation)

# ### A pretty movie
#
# We load the saved output and make a pretty movie of the simulation. First we plot a snapshot:
using CairoMakie

u = FieldTimeSeries("global_surface_fields_with_bgc.jld2", "u"; backend = OnDisk())
v = FieldTimeSeries("global_surface_fields_with_bgc.jld2", "v"; backend = OnDisk())
T = FieldTimeSeries("global_surface_fields_with_bgc.jld2", "T"; backend = OnDisk())
e = FieldTimeSeries("global_surface_fields_with_bgc.jld2", "e"; backend = OnDisk())

C = FieldTimeSeries("global_surface_fields_with_bgc.jld2", "DIC"; backend = OnDisk())
p = FieldTimeSeries("global_surface_fields_with_bgc.jld2", "ocean_pCO₂"; backend = OnDisk())
f = FieldTimeSeries("global_surface_fields_with_bgc.jld2", "CO₂_flux"; backend = OnDisk())

times = u.times
Nt = length(times)

n = Observable(Nt)

# We create a land mask and use it to fill land points with `NaN`s.
land = interior(T.grid.immersed_boundary.bottom_height) .≥ 0

Tn = @lift begin
    Tn = interior(T[$n])
    Tn[land] .= NaN
    view(Tn, :, :, 1)
end

en = @lift begin
    en = interior(e[$n])
    en[land] .= NaN
    view(en, :, :, 1)
end

Cn = @lift begin
    Cn = interior(C[$n])
    Cn[land] .= NaN
    view(Cn, :, :, 1)
end

pn = @lift begin
    pn = interior(p[$n])
    pn[land] .= NaN
    view(pn, :, :, 1)
end

fn = @lift begin
    fn = interior(f[$n])
    fn[land] .= NaN
    view(fn, :, :, 1)
end

# We compute the surface speed.
un = Field{Face, Center, Nothing}(u.grid)
vn = Field{Center, Face, Nothing}(v.grid)
s = Field(sqrt(un^2 + vn^2))

sn = @lift begin
    parent(un) .= parent(u[$n])
    parent(vn) .= parent(v[$n])
    compute!(s)
    sn = interior(s)
    sn[land] .= NaN
    view(sn, :, :, 1)
end

# Finally, we plot a snapshot of the surface speed, temperature, and the turbulent
# eddy kinetic energy from the CATKE vertical mixing parameterization.
fig = Figure(size = (800, 1200))

axs = Axis(fig[1, 1], xlabel="Longitude (deg)", ylabel="Latitude (deg)")
axT = Axis(fig[2, 1], xlabel="Longitude (deg)", ylabel="Latitude (deg)")
axe = Axis(fig[3, 1], xlabel="Longitude (deg)", ylabel="Latitude (deg)")

hm = heatmap!(axs, sn, colorrange = (0, 0.5), colormap = :deep, nan_color=:lightgray)
Colorbar(fig[1, 2], hm, label = "Surface speed (m s⁻¹)")

hm = heatmap!(axT, Tn, colorrange = (-1, 30), colormap = :magma, nan_color=:lightgray)
Colorbar(fig[2, 2], hm, label = "Surface Temperature (ᵒC)")

hm = heatmap!(axe, en, colorrange = (0, 1e-3), colormap = :solar, nan_color=:lightgray)
Colorbar(fig[3, 2], hm, label = "Turbulent Kinetic Energy (m² s⁻²)")

save("global_snapshot.png", fig)
nothing #hide

# ![](global_snapshot.png)

# And now a movie:

record(fig, "one_degree_global_ocean_surface.mp4", 1:Nt, framerate = 8) do nn
    n[] = nn
end
nothing #hide

fig = Figure(size = (800, 1200))

axC = Axis(fig[1, 1], xlabel="Longitude (deg)", ylabel="Latitude (deg)")
axp = Axis(fig[2, 1], xlabel="Longitude (deg)", ylabel="Latitude (deg)")
axf = Axis(fig[3, 1], xlabel="Longitude (deg)", ylabel="Latitude (deg)")

hm = heatmap!(axC, Cn, colorrange = (1.8, 2.2), colormap = :viridis, nan_color=:lightgray)
Colorbar(fig[1, 2], hm, label = "Surface DIC concentration (mol m⁻³)")

hm = heatmap!(axp, pn, colorrange = (300e-6, 480e-6), colormap = Reverse(:RdYlBu), nan_color=:lightgray)
Colorbar(fig[2, 2], hm, label = "Surface ocean pCO₂ (atm)")

hm = heatmap!(axf, fn, colorrange = (-5e-9, 5e-9), colormap = :bluesreds, nan_color=:lightgray)
Colorbar(fig[3, 2], hm, label = "Sea-air CO₂ fluxes (x10⁻⁷ mol m⁻² s⁻¹)")
save("snapshot_bgc.png", fig)
nothing #hide

# ![](snapshot.png)

# And now a movie:

record(fig, "near_global_ocean_surface_bgc.mp4", 1:Nt, framerate = 8) do nn
    n[] = nn
end
nothing #hide

# ![](near_global_ocean_surface_bgc.mp4)


