using ClimaOcean
using ClimaOcean.ECCO: ECCO4Monthly, ECCO4DarwinMonthly, ECCOFieldTimeSeries, NearestNeighborInpainting, interpolate_to_grid
using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients
using ClimaOceanBiogeochemistry.CarbonSystemSolvers.UniversalRobustCarbonSolver: UniversalRobustCarbonSystem
using ClimaOceanBiogeochemistry.CarbonSystemSolvers: CarbonSystemParameters, CarbonSolverParameters, CarbonCoefficientParameters
using OrthogonalSphericalShellGrids
using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: Time
using CFTime
using Dates
using Printf
using CUDA: @allowscalar, device!

using Oceananigans.Grids: znode

include("calculate_air_sea_carbon_exchange.jl")

arch = CPU()

#####
##### Grid and Bathymetry
#####

Nx = 360
Ny = 180
Nz = 100

z_faces = exponential_z_faces(; Nz, depth=5000, h=34)

underlying_grid = TripolarGrid(arch;
                               size = (Nx, Ny, Nz),
                               z = z_faces,
                               first_pole_longitude = 70,
                               north_poles_latitude = 55)

bottom_height = regrid_bathymetry(underlying_grid;
                                  minimum_depth = 10,
                                  interpolation_passes = 75,
                                  major_basins = 2)

# Open Gibraltar strait 
# TODO: find a better way to do this
tampered_bottom_height = deepcopy(bottom_height)
view(tampered_bottom_height, 102:103, 124, 1) .= -400

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(tampered_bottom_height))

#####
##### Closures
#####

gm = Oceananigans.TurbulenceClosures.IsopycnalSkewSymmetricDiffusivity(κ_skew=1000, κ_symmetric=1000)
catke = ClimaOcean.OceanSimulations.default_ocean_closure()
viscous_closure = Oceananigans.TurbulenceClosures.HorizontalScalarDiffusivity(ν=2000)

closure = (gm, catke, viscous_closure)

#####
##### Restoring
#####

restoring_rate  = 1 / 2days
z_below_surface = @allowscalar znode(1, 1, grid.Nz, grid, Center(), Center(), Face())

mask = LinearlyTaperedPolarMask(southern=(-80, -70), northern=(70, 90), z=(z_below_surface, 0))

start = DateTimeProlepticGregorian(1993, 1, 1)
stop  = DateTimeProlepticGregorian(1993, 12, 1)
dates = range(start; stop, step=Month(1))
temperature = ECCOMetadata(:temperature; dates, version=ECCO4Monthly(), dir="./")
salinity    = ECCOMetadata(:salinity;    dates, version=ECCO4Monthly(), dir="./")

# inpainting = NearestNeighborInpainting(30) should be enough to fill the gaps near bathymetry
FT = ECCORestoring(:temperature, grid; dates=dates, mask=mask, rate=restoring_rate, inpainting=NearestNeighborInpainting(50))
FS = ECCORestoring(:salinity,    grid; dates=dates, mask=mask, rate=restoring_rate, inpainting=NearestNeighborInpainting(50))

#####
##### Ocean simulation
##### 

momentum_advection = VectorInvariant()
tracer_advection   = Centered(order=2)

# biogeochem
# store boundary condition for DIC
CO₂_flux   = Field{Center, Center, Nothing}(grid)

## These are filled in compute_CO₂_flux!
ocean_pCO₂ = Field{Center, Center, Nothing}(grid)
atmos_pCO₂ = Field{Center, Center, Nothing}(grid)
pH         = Field{Center, Center, Nothing}(grid)
set!(pH, 8.0)

## This is filled in set_par!
PAR       = Field{Center, Center, Nothing}(grid)

biogeochemistry = CarbonAlkalinityNutrients(; 
            grid,
            maximum_net_community_production_rate = 0.5/365.25days,
            incident_PAR = PAR,
)

carbon     = ECCOMetadata(:DIC; dates, version=ECCO4DarwinMonthly(), dir="./")
alkalinity = ECCOMetadata(:ALK; dates, version=ECCO4DarwinMonthly(), dir="./")
phosphate  = ECCOMetadata(:PO₄; dates, version=ECCO4DarwinMonthly(), dir="./")
nitrate    = ECCOMetadata(:NO₃; dates, version=ECCO4DarwinMonthly(), dir="./")
iron       = ECCOMetadata(:Fe;  dates, version=ECCO4DarwinMonthly(), dir="./")
dissolvedOganicPhosphate = ECCOMetadata(
                                :DOP; dates, version=ECCO4DarwinMonthly(), dir="./",
)
particulateOrganicPhosphate = ECCOMetadata(
                                :POP; dates, version=ECCO4DarwinMonthly(), dir="./",
)
# For pCO2 calculation
silicate = ECCOMetadata(:Siᵀ; dates, version=ECCO4DarwinMonthly(), dir="./")
FSiᵀ     = ECCOFieldTimeSeries(silicate, grid)

# Initialize carbon surface boundary condition (CO2 fluxes)
BC  = FieldBoundaryConditions(
    top = FluxBoundaryCondition(CO₂_flux) #+FluxBoundaryCondition(VirtualFlux)
    )
boundary_conditions = (DIC=BC,) #, ALK=FA)

# Initialatize the carbon and alkalinity FW forcing (related to salinity restoring)
FC = ECCORestoring(:salinity,    grid; dates=dates, mask=mask, rate=2.2*(restoring_rate/34.5), inpainting=NearestNeighborInpainting(50))
FA = ECCORestoring(:salinity,    grid; dates=dates, mask=mask, rate=2.4*(restoring_rate/34.5), inpainting=NearestNeighborInpainting(50))
#FF = FieldBoundaryConditions(
#    top = FluxBoundaryCondition(DustDeposition)
#    )

forcing = (T=FT, S=FS, DIC=FC, ALK=FA) # , FE=FF)

# Should we add a side drag since this is at a coarser resolution?
ocean = ocean_simulation(grid; momentum_advection, tracer_advection,
                         closure, forcing, boundary_conditions, biogeochemistry,
                         tracers = (:T, :S, :e, :DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe))

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

#####
##### Atmospheric forcing
#####

radiation  = Radiation(arch)
atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(20))

# Set the PAR to 40% of the downwelling shortwave radiation
set!(PAR, 0.4 * [interpolate_to_grid(
        atmosphere.downwelling_radiation.shortwave, 
        i, j, 1, 
        atmosphere.downwelling_radiation.shortwave.grid, 
        grid, 
        Oceananigans.Units.Time(0),
        ) for i in 1:Nx, j in 1:Ny],
)

@kernel function interpolate_atm2oce!(
        grid::ImmersedBoundaryGrid, 
        ocean_forcing::Field, 
        atm_forcing::FieldTimeSeries, 
        mytime::Oceananigans.Units.Time, 
        scale_frac::Real = 1,
        )   
    i, j = @index(Global, NTuple)
    k = size(grid, 3)
    inactive = inactive_cell(i, j, k, grid)
        
    while !inactive
        @inbounds begin    
            ocean_forcing[i, j, 1] = scale_frac * interpolate_to_grid(
                atm_forcing, 
                i, j, 1, 
                atm_forcing.grid, 
                grid, 
                mytime,
                )
        end
    end
end

@inline function update_par!(simulation::Simulation, PAR_frac::Real = 0.4 )
    kernel_args = (
        simulation.model.ocean.model.grid,
        simulation.model.ocean.model.biogeochemistry.incident_PAR,
        simulation.model.atmosphere.downwelling_radiation.shortwave,
        Time(simulation.model.clock),
        PAR_frac,
    )

    launch!(arch,
        simulation.model.ocean.model.grid,
        :xy,
        interpolate_atm2oce!,
        kernel_args...,
    )
    return nothing
end

#####
##### Coupled simulation
#####

coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation) 
simulation    = Simulation(coupled_model; Δt=15minutes, stop_time=30minutes) #2*365days)

#####
##### Run it!
##### 

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

    @info @sprintf("Time: %s, n: %d, Δt: %s, max|u|: (%.2e, %.2e, %.2e) m s⁻¹, extrema(T): (%.2f, %.2f) ᵒC, wall time: %s \n",
                   prettytime(sim), iteration(sim), prettytime(sim.Δt),
                   umax..., Tmax, Tmin, prettytime(step_time))

     wall_time[] = time_ns()

     return nothing
end

add_callback!(simulation, progress, IterationInterval(10))
add_callback!(simulation, update_par!)
add_callback!(simulation, calculate_air_sea_carbon_exchange!)
run!(simulation)
