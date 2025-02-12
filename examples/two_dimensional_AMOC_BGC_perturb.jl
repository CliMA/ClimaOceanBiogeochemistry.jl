# Pick up from a baseline model, run 1-year perturbation
using GLMakie
using Printf
using Statistics

using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients
using Oceananigans
using Oceananigans.Models.HydrostaticFreeSurfaceModels:
                    HydrostaticFreeSurfaceModel,
                    PrescribedVelocityFields
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization
using Oceananigans.Units
using Oceananigans.Fields: ZeroField, CenterField, FunctionField
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans: TendencyCallsite

Ny = 500 
Nz = 200
Ly = 15000kilometers   # m
Lz = 4000           # m

arch = CPU()
# We use a two-dimensional grid, with a `Flat` `y`-direction:
grid = RectilinearGrid(arch,
                       size = (Ny, Nz),
                       y = (0, Ly),
                       z = (-Lz, 0),
                       topology=(Flat, Bounded, Bounded))

############################## Physical circulation ############################## 
# Set streamfunction
deltaN = 500kilometers   # North downwelling width
deltaS = 3000kilometers  # South upwelling width
deltaZ = 1000     # Vertical asymmetry
Ψᵢ(y, z)  = - 7 * ((1 - exp(-y / deltaS)) * (1 - exp(-(Ly - y) / deltaN)) * 
            sinpi(z / Lz) * exp(z/deltaZ))

Ψ = Field{Center, Face, Face}(grid)
set!(Ψ, Ψᵢ)
fill_halo_regions!(Ψ, arch)

# Set velocity field from streamfunction
v = YFaceField(grid)
w = ZFaceField(grid)
v .= - ∂z(Ψ)
w .= + ∂y(Ψ)
fill_halo_regions!(v, arch)
fill_halo_regions!(w, arch)

# (I have to specify u to allow CheckPointer)
u = XFaceField(grid) 
fill_halo_regions!(u, arch)

############################# Diffusivity ############################# 

# Seasonal: 100+50*sinpi(2*(t/day)/365.25)
# Single: (t % 10days == 0 ? 200 : 0)
# Initialize the perturbation after Day 1 (not Day 0): 365.25*2000 + 1
kz(y,z,t) = 1e-4 + 5e-3 * (tanh((z+(100+50*sinpi(2*(t/day-(365.25*2000 + 1))/(365/12))))/20)+1) + 1e-2 * exp(-(z+4000)/50)
tracer_vertical_closure = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), 
                                        κ=(DIC=kz,ALK=kz,PO₄=kz,NO₃=kz,DOP=kz,POP=kz,Fe=0))
tracer_horizontal_closure = HorizontalScalarDiffusivity(
                                    κ=(DIC=1e3,ALK=1e3,PO₄=1e3,NO₃=1e3,DOP=1e3,POP=1e3,Fe=0))

# Fe field & forcing
# @inline Feᵢ(y,z) = (0.02+(-z/2000)) /1e6
# Fe_input(y, z, t) = (z > -30 && y < 100kilometers ? 1e-12 : 0.0)
# F_forcing = Forcing(Fe_input)

maximum_net_community_production_rate = CenterField(grid) 
maxNCP(y,z) = (2e-5 + sinpi(y/Ly) * 8e-5)/day
set!(maximum_net_community_production_rate, maxNCP)   
fill_halo_regions!(maximum_net_community_production_rate, arch)

############################### Set PAR(y,t) ############################### 
#=
function seasonal_PAR(y, z, t)
    # Constants
    day_in_year = 365.25
    axial_tilt = 23.5  # Earth's axial tilt in degrees
    latitude = (y/Ly - 0.5) * 180.0  # Convert y to latitude in degrees
    
    # Solar declination (varies over the year due to Earth's tilt)
    solar_declination = axial_tilt * sinpi(2 * t / day_in_year)
    
    # Solar angle of incidence (latitude and solar declination combined)
    solar_angle = latitude - solar_declination
    
    # PAR depends on the cosine of the solar angle 
    angle_effect = cosd(solar_angle)
    angle_effect = max(0.0, angle_effect)  # Ensure no negative values
    
    # Scale PAR by angle effect
    PAR = 700 * angle_effect
    return PAR
end

clock = Clock{Float64}(time=0)
incident_PAR = FunctionField{Nothing, Center, Center}(seasonal_PAR, grid; clock)
=#
# Set PAR as a function of latitude
incident_PAR = CenterField(grid) 
surface_PAR(y,z) = 700 * sinpi(y/Ly) 
set!(incident_PAR, surface_PAR)   
fill_halo_regions!(incident_PAR, arch)

############################## Model 2: Perturbation ############################## 
model2 = HydrostaticFreeSurfaceModel(grid = grid,
                                     #clock,
                                     biogeochemistry = CarbonAlkalinityNutrients(; grid,
                                                                                   maximum_net_community_production_rate  = maximum_net_community_production_rate,
                                                                                   incident_PAR = incident_PAR,
                                                                                   particulate_organic_phosphorus_sedremin_timescale = 0.5 / day),
                                                                                   #iron_scavenging_rate = 0),
                                     velocities = PrescribedVelocityFields(; u, v, w),
                                     tracers = (:DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe),
                                     # forcing=(; Fe=F_forcing),
                                     coriolis = nothing,
                                     buoyancy = nothing,
                                     closure = (tracer_vertical_closure, tracer_horizontal_closure))

set!(model2, DIC=2.1, ALK=2.35, NO₃=2.4e-2, PO₄=1.6e-3, DOP=0, POP=0, Fe = 6e-7) # mol PO₄ m⁻³

spinup_time = 365.25*2000days
perturbation_time =366days
simulation2 = Simulation(model2; Δt = 1days, stop_time=spinup_time+perturbation_time) 

# Define a callback to zero out Fe tendency
function modify_tendency!(model2)
        model2.timestepper.Gⁿ.Fe .= 0
        # model2.timestepper.Gⁿ.Fe[1,1:50,100] .= 1e-11
        return nothing
end                                        
simulation2.callbacks[:modify_Fe] = Callback(modify_tendency!, 
                                            callsite = TendencyCallsite())

outputs = (v = model2.velocities.v,
            w = model2.velocities.w,
            PO₄= model2.tracers.PO₄,
            DOP = model2.tracers.DOP,
            POP = model2.tracers.POP,
            NO₃ = model2.tracers.NO₃,
            Fe = model2.tracers.Fe,
            NCP = model2.biogeochemistry.NCP,
            Premin = model2.biogeochemistry.Premin,
            Dremin = model2.biogeochemistry.Dremin)

simulation2.output_writers[:simple_output] =
            JLD2OutputWriter(model2, outputs, 
                            schedule = TimeInterval(1days), 
                            filename = "AMOC115_12perY",
                            overwrite_existing = true)   

simulation2.output_writers[:checkpointer] = Checkpointer(model2,
            schedule = TimeInterval(perturbation_time),
            prefix = "AMOC115_checkpoint",
            overwrite_existing = false)
            
run!(simulation2, pickup = true)

#################################### Video of all tracers ####################################

filepath = simulation2.output_writers[:simple_output].filepath
# filepath = "./AMOC115_seasonal.jld2"

PO4_timeseries = FieldTimeSeries(filepath, "PO₄")
times = PO4_timeseries.times
xw, yw, zw = nodes(PO4_timeseries)
POP_timeseries = FieldTimeSeries(filepath, "POP")
DOP_timeseries = FieldTimeSeries(filepath, "DOP")
Fe_timeseries = FieldTimeSeries(filepath, "Fe")
NO3_timeseries = FieldTimeSeries(filepath, "NO₃")

NCP_timeseries = FieldTimeSeries(filepath, "NCP")
Premin_timeseries = FieldTimeSeries(filepath, "Premin")
Dremin_timeseries = FieldTimeSeries(filepath, "Dremin")

#=
n = Observable(1)
# title = @lift @sprintf("t = Year %d0", times[$n] / 3652.5days) 
title = @lift @sprintf("t = Day %d0", times[$n] / 10days) 

# convert unit from mol/m³ to μM: 1e3*interior(...)
PO4ₙ = @lift 1e3*interior(PO4_timeseries[$n], 1, :, :)
avg_PO4ₙ = @lift mean(1e3*interior(PO4_timeseries[$n], 1, :, :), dims=1) 

POPₙ = @lift 1e3*interior(POP_timeseries[$n], 1, :, :)
avg_POPₙ = @lift mean(1e3*interior(POP_timeseries[$n], 1, :, :), dims=1) 

NO3ₙ = @lift 1e3*interior(NO3_timeseries[$n], 1, :, :)
Feₙ = @lift 1e6*interior(Fe_timeseries[$n], 1, :, :)

fig = Figure(size = (1500, 1500))

ax_s = Axis(fig[2, 1]; xlabel = "y (km)", ylabel = "z (m)",title = "Stream function (m² s⁻¹)", aspect = 1)
hm_s = heatmap!(ax_s, yw./1e3, zw, Ψ; colormap = :viridis) 
Colorbar(fig[2, 2], hm_s; flipaxis = false)
contour!(ax_s, yw./1e3, grid.zᵃᵃᶜ[0:100], Ψ, levels = 10, color = :black)

ax_v = Axis(fig[2, 3]; xlabel = "y (km)", ylabel = "z (m)",title = "v (t) (cm s⁻¹)", aspect = 1)
hm_v = heatmap!(ax_v, yw./1e3, zw, vₙ; colorrange = (-3,3), colormap = :balance) 
Colorbar(fig[2, 4], hm_v; flipaxis = false)

ax_w = Axis(fig[2, 5]; xlabel = "y (km)", ylabel = "z (m)", title = "w (t) (cm s⁻¹)", aspect = 1)
hm_w = heatmap!(ax_w, yw./1e3, zw, wₙ; colorrange = (-8e-4, 8e-4), colormap = :balance) 
Colorbar(fig[2, 6], hm_w; flipaxis = false)

ax_PO4 = Axis(fig[3, 1]; xlabel = "y (km)", ylabel = "z (m)", title = "[PO₄] (μM)", aspect = 1)
hm_PO4 = heatmap!(ax_PO4, yw/1e3, zw, PO4ₙ; colorrange = (0,2.5),colormap = :rainbow1) 
Colorbar(fig[3, 2], hm_PO4; flipaxis = false)
contour!(ax_PO4, yw/1e3, zw[1:100], PO4ₙ, levels = 5, color = :black)

ax_NO3 = Axis(fig[3, 3]; xlabel = "y (km)", ylabel = "z (m)", title = "[NO₃] (μM)", aspect = 1)
hm_NO3 = heatmap!(ax_NO3, yw/1e3, zw, NO3ₙ; colorrange = (0,35),colormap = :rainbow1) 
Colorbar(fig[3, 4], hm_NO3; flipaxis = false)
contour!(ax_NO3, yw/1e3, zw[1:100], NO3ₙ, levels = 5, color = :black)

ax_Fe = Axis(fig[3, 5]; xlabel = "y (km)", ylabel = "z (m)", title = "[Fe] (nM)", aspect = 1)
hm_Fe = heatmap!(ax_Fe, yw/1e3, zw, Feₙ; colorrange = (0,1),colormap = :rainbow1) 
Colorbar(fig[3, 6], hm_Fe; flipaxis = false)
contour!(ax_Fe, yw/1e3, zw[1:100], Feₙ, levels = 5, color = :black)

ax_POP = Axis(fig[4, 3]; xlabel = "y (km)", ylabel = "z (m)", title = "[POP] (μM)", aspect = 1)
hm_POP = heatmap!(ax_POP, yw/1e3, zw, POPₙ; colorrange = (0,0.015),colormap = :rainbow1) 
Colorbar(fig[4, 4], hm_POP; flipaxis = false)

ax_avg_PO4 = Axis(fig[4, 1:2]; xlabel = "[PO₄] (μM)", ylabel = "z (m)", title = "Average [PO₄] (μM)", yaxisposition = :right)
xlims!(ax_avg_PO4, 0, 2)
PO4_prof = lines!(ax_avg_PO4, avg_PO4ₙ[][1, :], zw[1:100])

ax_avg_POP = Axis(fig[4, 5:6]; xlabel = "[POP] (μM)", ylabel = "z (m)", title = "Average [POP] (μM)",yaxisposition = :right)
xlims!(ax_avg_POP, 0, 0.015)
POP_prof = lines!(ax_avg_POP, avg_POPₙ[][1, :], zw[1:100])

fig[1, 1:6] = Label(fig, title, tellwidth=false)

# And, finally, we record a movie.
frames = 1:length(times)
record(fig, "AMOC28_perturb.mp4", frames, framerate=40) do i
    n[] = i
    PO4_prof[1] = avg_PO4ₙ[][1, :]
    POP_prof[1] = avg_POPₙ[][1, :]
end
nothing #hide
=#
############################## Compare concentration variatons ##############################

# PO4_last = 1e3*interior(PO4_timeseries[end], 1, :, :) 
# POP_last = 1e3*interior(POP_timeseries[end], 1, :, :) 
# DOP_last = 1e3*interior(DOP_timeseries[end], 1, :, :) 
# sum(PO4_last.+POP_last.+DOP_last)

n = Observable(1)
title = @lift @sprintf("t = Day %d", ((times[$n] - 365.25*2000days) / 1days)-1) 

PO4_init = 1e3*interior(PO4_timeseries[1], 1, :, :) 
POP_init = 1e3*interior(POP_timeseries[1], 1, :, :) 
DOP_init = 1e3*interior(DOP_timeseries[1], 1, :, :)

diff_PO4 = @lift 1e3*(interior(PO4_timeseries[$n], 1, :, :) .- interior(PO4_timeseries[1], 1, :, :))
diff_POP = @lift 1e3*(interior(POP_timeseries[$n], 1, :, :) .- interior(POP_timeseries[1], 1, :, :))
diff_DOP = @lift 1e3*(interior(DOP_timeseries[$n], 1, :, :) .- interior(DOP_timeseries[1], 1, :, :))

# flux = (10m/d)*[POP]
avg_POP_flux_init = mean(1e4*interior(POP_timeseries[1], 1, :, :) , dims=1) 
avg_POP_fluxₙ = @lift mean(1e4*interior(POP_timeseries[$n], 1, :, :), dims=1) 

POP_flux_init = 1e4*interior(POP_timeseries[1], 1, :, :) 
POP_fluxₙ = @lift 1e4*interior(POP_timeseries[$n], 1, :, :)

z₀ = log(0.01)*25 
Martin_factor = ((zw[190] .+ z₀) ./ (zw .+ z₀)).^0.84

# 1D avergae profile
avg_Martin_fluxₙ = @lift($avg_POP_fluxₙ[190] .* Martin_factor)

# 2D heatmap
Martin_fluxₙ = Observable(similar(POP_fluxₙ[]))  

Martin_fluxₙ = @lift begin
    matrix = Observable(similar($POP_fluxₙ))
    for i in 1:200
        matrix[][:, i] .= $POP_fluxₙ[:, 190] .* $Martin_factor[i]
    end
    return matrix
end
diff_Martin = @lift ($POP_fluxₙ .- $Martin_fluxₙ[])

# One station
# up_PO4ₙ = @lift 1e3*interior(PO4_timeseries[$n], 1, 250, :)
# up_DOPₙ = @lift 1e3*interior(DOP_timeseries[$n], 1, 250, :)
# up_POPₙ = @lift 1e3*interior(POP_timeseries[$n], 1, 250, :)

# init_up_PO4ₙ = 1e3*interior(PO4_timeseries[1], 1, 250, :)
# init_up_POPₙ = 1e3*interior(POP_timeseries[1], 1, 250, :)
# init_up_DOPₙ = 1e3*interior(DOP_timeseries[1], 1, 250, :)


t= collect(0:1days:perturbation_time)
kz_profiles = Matrix{Float64}(undef, length(zw), length(t))
for (i, t_val) in enumerate(t)
    kz_profiles[:, i] = [kz(0, z_val, t_val) for z_val in zw]
end

fig_compare = Figure(size=(1000, 1000))

ax_kz = Axis(fig_compare[2, 1:2], xlabel = "kz (m² s⁻¹)", ylabel = "z (m)", title = "Temporal Variation of kz")
kz_prof = lines!(ax_kz, kz_profiles[:, 1], zw)
ylims!(ax_kz, -1000, 0)

# Column 1: tracers final - initial
ax_PO4 = Axis(fig_compare[2, 3]; xlabel = "y (km)", ylabel = "z (m)", title = "Δ[PO₄] (μM)", aspect = 1)
hm_PO4 = heatmap!(ax_PO4, yw/1e3, zw, diff_PO4; colorrange = (-0.2,0.2),colormap = :balance) 
Colorbar(fig_compare[2, 4], hm_PO4; flipaxis = false)
ylims!(ax_PO4, -1000, 0)

ax_POP = Axis(fig_compare[3, 3]; xlabel = "y (km)", ylabel = "z (m)", title = "Δ[POP] (μM)", aspect = 1)
hm_POP = heatmap!(ax_POP, yw/1e3, zw, diff_POP; colorrange = (-0.002,0.002),colormap = :balance) 
Colorbar(fig_compare[3, 4], hm_POP; flipaxis = false)
ylims!(ax_POP, -1000, 0)

ax_DOP = Axis(fig_compare[4, 3]; xlabel = "y (km)", ylabel = "z (m)", title = "Δ[DOP] (μM)", aspect = 1)
hm_DOP = heatmap!(ax_DOP, yw/1e3, zw, diff_DOP; colorrange = (-0.07,0.07),colormap = :balance) 
Colorbar(fig_compare[4, 4], hm_DOP; flipaxis = false)
ylims!(ax_DOP, -1000, 0)

# POP flux 
ax_avg_POP_flux = Axis(fig_compare[3, 1:2]; xlabel = "POP flux (mmol m⁻² d⁻¹)", ylabel = "z (m)", title = "Average POP flux")
lines!(ax_avg_POP_flux, vec(avg_POP_flux_init), zw, color = :grey, linewidth = 1, label = "Initial")
POP_flux_prof = lines!(ax_avg_POP_flux, avg_POP_fluxₙ[][1,:], zw,  linewidth = 2, label = "Dynamic")
Martin_flux_prof = lines!(ax_avg_POP_flux, avg_Martin_fluxₙ[], zw, color = :red,  linewidth = 2, linestyle = :dash, label = "Martin")
ylims!(ax_avg_POP_flux, -1000, 0)
xlims!(ax_avg_POP_flux, 0, 0.08)
axislegend(ax_avg_POP_flux, position = :rb)

# diff POP flux - Martin
ax_diffMartin = Axis(fig_compare[4, 1]; xlabel = "y (km)", ylabel = "z (m)", title = "Δ(POP flux): True - Martin", aspect = 1)
hm_diffMartin = heatmap!(ax_diffMartin, yw/1e3, zw, diff_Martin; colorrange = (-7e-3,7e-3),colormap = :balance) 
Colorbar(fig_compare[4, 2], hm_diffMartin; flipaxis = false)
ylims!(ax_diffMartin, -1000, -200)

fig_compare[1, 1:4] = Label(fig_compare, title, tellwidth=false)

# And, finally, we record a movie.
frames = 1:length(times)
record(fig_compare, "AMOC115_12perY_conc.mp4", frames, framerate=50) do i
    n[] = i
    kz_prof[1]= kz_profiles[:,i]
    POP_flux_prof[1] = avg_POP_fluxₙ[][1,:]
    Martin_flux_prof[1] = avg_Martin_fluxₙ[]
end
nothing #hide

######################## Compare rate variatons ####################

NCP_init = Observable(1days*1e3*interior(NCP_timeseries[1], 1, :, :) )
Premin_init = Observable(1days*1e3*interior(Premin_timeseries[1], 1, :, :))
Dremin_init = Observable(1days*1e3*interior(Dremin_timeseries[1], 1, :, :) )

n = Observable(1)
title = @lift @sprintf("t = Day %d", ((times[$n] - 365.25*2000days) / 1days)-1)  

# Difference
diff_NCP = @lift (1days*1e3).*(interior(NCP_timeseries[$n], 1, :, :) .- interior(NCP_timeseries[1], 1, :, :))
diff_Premin = @lift (1days*1e3).*(interior(Premin_timeseries[$n], 1, :, :) .- interior(Premin_timeseries[1], 1, :, :))
diff_Dremin = @lift (1days*1e3).*(interior(Dremin_timeseries[$n], 1, :, :) .- interior(Dremin_timeseries[1], 1, :, :))

# Anomaly
anomaly_NCP = @lift (($diff_NCP ./ $NCP_init) .* 100)
anomaly_Premin = @lift (($diff_Premin ./ $Premin_init) .* 100)
anomaly_Dremin = @lift (($diff_Dremin ./ $Dremin_init) .* 100)

# Single station
up_NCPₙ = @lift 1days*1e3*interior(NCP_timeseries[$n], 1, 250, :)
up_Preminₙ = @lift 1days*1e3*interior(Premin_timeseries[$n], 1, 250, :)
up_Dreminₙ = @lift 1days*1e3*interior(Dremin_timeseries[$n], 1, 250, :)

init_up_NCPₙ = 1days*1e3*interior(NCP_timeseries[1], 1, 250, :)
init_up_Preminₙ = 1days*1e3*interior(Premin_timeseries[1], 1, 250, :)
init_up_Dreminₙ = 1days*1e3*interior(Dremin_timeseries[1], 1, 250, :)

t= collect(0:1days:perturbation_time)
kz_profiles = Matrix{Float64}(undef, length(zw), length(t))
for (i, t_val) in enumerate(t)
    kz_profiles[:, i] = [kz(0, z_val, t_val) for z_val in zw]
end

### Plot ###
fig_rates = Figure(size=(1500, 1200))

ax_kz = Axis(fig_rates[2, 1], xlabel = "kz (m² s⁻¹)", ylabel = "z (m)", title = "Temporal Variation of kz")
kz_prof = lines!(ax_kz, kz_profiles[:, 1], zw)
ylims!(ax_kz, -1000, 0)

# Column 2: rates final - initial
ax_NCP = Axis(fig_rates[2, 2]; xlabel = "y (km)", ylabel = "z (m)", title = "Δ(NCP) (mmol P m⁻³ d⁻¹)", aspect = 1)
hm_NCP = heatmap!(ax_NCP, yw/1e3, zw, diff_NCP; colorrange = (-1e-3,1e-3),colormap = :balance) 
Colorbar(fig_rates[2, 3], hm_NCP; flipaxis = false)
ylims!(ax_NCP, -1000, 0)

ax_Premin = Axis(fig_rates[3, 2]; xlabel = "y (km)", ylabel = "z (m)", title = "Δ(POPremin) (mmol P m⁻³ d⁻¹)", aspect = 1)
hm_Premin = heatmap!(ax_Premin, yw/1e3, zw, diff_Premin; colorrange = (-1e-4,1e-4),colormap = :balance) 
Colorbar(fig_rates[3, 3], hm_Premin; flipaxis = false)
ylims!(ax_Premin, -1000, 0)

ax_Dremin = Axis(fig_rates[4, 2]; xlabel = "y (km)", ylabel = "z (m)", title = "Δ(DOPremin) (mmol P m⁻³ d⁻¹)", aspect = 1)
hm_Dremin = heatmap!(ax_Dremin, yw/1e3, zw, diff_Dremin; colorrange = (-5e-4,5e-4),colormap = :balance) 
Colorbar(fig_rates[4, 3], hm_Dremin; flipaxis = false)
ylims!(ax_Dremin, -1000, 0)

# Column 3: rates anomaly
ax_anomalyNCP = Axis(fig_rates[2, 4]; xlabel = "y (km)", ylabel = "z (m)", title = "NCP anomaly (% diff/init)", aspect = 1)
hm_anomalyNCP = heatmap!(ax_anomalyNCP, yw/1e3, zw, anomaly_NCP; colorrange = (-100,100), colormap = :balance) 
Colorbar(fig_rates[2, 5], hm_anomalyNCP; flipaxis = false)
ylims!(ax_anomalyNCP, -1000, 0)

ax_anomalyPremin = Axis(fig_rates[3, 4]; xlabel = "y (km)", ylabel = "z (m)", title = "POP remin anomaly (% diff/init)", aspect = 1)
hm_anomalyPremin = heatmap!(ax_anomalyPremin, yw/1e3, zw, anomaly_Premin; colorrange = (-100,100), colormap = :balance) 
Colorbar(fig_rates[3, 5], hm_anomalyPremin; flipaxis = false)
ylims!(ax_anomalyPremin, -1000, 0)

ax_anomalyDremin = Axis(fig_rates[4, 4]; xlabel = "y (km)", ylabel = "z (m)", title = "DOP remin anomaly (% diff/init)", aspect = 1)
hm_anomalyDremin = heatmap!(ax_anomalyDremin, yw/1e3, zw, anomaly_Dremin; colorrange = (-100,100), colormap = :balance) 
Colorbar(fig_rates[4, 5], hm_anomalyDremin; flipaxis = false)
ylims!(ax_anomalyDremin, -1000, 0)

# Column 4: rates: single profile 
up_NCP = Axis(fig_rates[2, 6]; title = "NCP (center station)", xlabel = "rate (mmol P m⁻³ d⁻¹)", ylabel = "z (m)")
lines!(up_NCP, init_up_NCPₙ, zw, color = :red, linestyle = :dash, label = "Initial")
NCP_up = lines!(up_NCP, up_NCPₙ[], zw, label = "Dynamic")
axislegend(up_NCP, position = :rb)
xlims!(up_NCP, 0, 5e-3)
ylims!(up_NCP, -1000, 0)

up_Premin = Axis(fig_rates[3, 6]; title = "POP remin (center station)", xlabel = "rate (mmol P m⁻³ d⁻¹)", ylabel = "z (m)")
lines!(up_Premin, init_up_Preminₙ, zw, color = :red, linestyle = :dash, label = "Initial")
Premin_up = lines!(up_Premin, up_Preminₙ[], zw, label = "Dynamic")
axislegend(up_Premin, position = :rb)
xlims!(up_Premin, 0, 5e-4)
ylims!(up_Premin, -1000, 0)

up_Dremin = Axis(fig_rates[4, 6]; title = "DOP remin (center station)", xlabel = "rate (mmol P m⁻³ d⁻¹)", ylabel = "z (m)")
lines!(up_Dremin, init_up_Dreminₙ, zw, color = :red, linestyle = :dash, label = "Initial")
Dremin_up = lines!(up_Dremin, up_Dreminₙ[], zw, label = "Dynamic")
axislegend(up_Dremin, position = :rb)
xlims!(up_Dremin, 0, 1.5e-3)
ylims!(up_Dremin, -1000, 0)

fig_rates[1, 1:6] = Label(fig_rates, title, tellwidth=false)

# And, finally, we record a movie.
frames = 1:length(times)
record(fig_rates, "AMOC115_12perY_rate.mp4", frames, framerate=50) do i
    n[] = i
    kz_prof[1]= kz_profiles[:,i]
    NCP_up[1] = up_NCPₙ[]
    Premin_up[1] = up_Preminₙ[]
    Dremin_up[1] = up_Dreminₙ[]
end
nothing #hide
#
######################## Fluxes vs Final: fit to Martin curve #############################
#=
end_down_POPₙ = 1e3*interior(POP_timeseries[end], 1, 300, :)
wₛₙₖ = 10 # m/day
z₀ = log(0.01)*25 

# initial and end flux: model
init_POP_flux = wₛₙₖ .* init_down_POPₙ
end_POP_flux = wₛₙₖ .* end_down_POPₙ

# initial and end flux: Martin
ref_grid_index = 191
init_Martin_flux = init_POP_flux[ref_grid_index]*((zw[ref_grid_index]+z₀)./(zw.+z₀)).^0.84
# init_Martin_flux[ref_grid_index+3:end] .= NaN

end_Martin_flux = end_POP_flux[ref_grid_index]*((zw[ref_grid_index]+z₀)./(zw.+z₀)).^0.84
# end_Martin_flux[ref_grid_index+3:end] .= NaN

fig_comp = Figure(size = (500, 500))

ax_flux = Axis(fig_comp[1, 1]; xlabel = "POP flux (mmol P m⁻² d⁻¹)", ylabel = "z (m)", title = "POP flux")
xlims!(ax_flux, 0, 0.1)
lines!(ax_flux, init_POP_flux[2:200], zw[2:200], linewidth = 3, color = :red, label = "Model (initial)")
lines!(ax_flux, init_Martin_flux, zw, linewidth = 4, color = :black, linestyle=:dash, label = "Martin (initial)")
lines!(ax_flux, end_POP_flux[2:200], zw[2:200], linewidth = 3, color = :blue, label = "Model (final)")
lines!(ax_flux, end_Martin_flux, zw, linewidth = 4, color = :grey, linestyle=:dash, label = "Martin (final)")
axislegend(ax_flux, position = :rb)
display(fig_comp)
=#

################################# Plot fluxes ################################# 
#=
time_after = 365
 
POP_t_up = zeros(Nz,time_after) # zw,n
POP_t_up_diff = zeros(Nz,time_after) 

# Initial peak flux: POP[97] ~80 m
# Final peak flux: POP[94] ~140 m
F100 = zeros(1,time_after) 
F140 = zeros(1,time_after) 
F200 = zeros(1,time_after) 
F500 = zeros(1,time_after) 
F1000 = zeros(1,time_after) 
F2000 = zeros(1,time_after) 

for i in 1:time_after
    POP_t_up[:, i] = 1e3*interior(POP_timeseries[i], 1, 10, :)
    POP_t_up_diff[:, i] = 1e3*(interior(POP_timeseries[i], 1, 10, :).-interior(POP_timeseries[1], 1, 10, :))
    
    F100[1,i] = 10 .* POP_t_up[96,i] 
    F140[1,i] = 10 .* POP_t_up[94,i] 
    F200[1,i] = 10 .* POP_t_up[91,i] #/ POP_t_up[93,i] 
    F500[1,i] = 10 .* POP_t_up[76,i] 
    F1000[1,i] = 10 .* POP_t_up[51,i] 
    F2000[1,i] = 10 .* POP_t_up[1,i] 
end
z₀ = -log(0.01)*25 
F200_martin = ones(1,time_after) .* ((200+z₀)/(100+z₀))^-0.84
F500_martin = ones(1,time_after) .* ((500+z₀)/(100+z₀))^-0.84
F1000_martin = ones(1,time_after) .* ((1000+z₀)/(100+z₀))^-0.84
F2000_martin = ones(1,time_after) .* ((2000+z₀)/(100+z₀))^-0.84

fig_flux = Figure(size=(800, 600))

ax_F1000 = Axis(fig_flux[1, 1]; xlabel = "time (days)", ylabel = "F(z) (mmol m⁻² d⁻¹)", yaxisposition = :left)
lines!(ax_F1000, 1:time_after, F140[1,:], color = :orange, label = "140 m (final peak)")
lines!(ax_F1000, 1:time_after, F200[1,:], color = :black, label = "200 m")
lines!(ax_F1000, 1:time_after, F500[1,:], color = :green, label = "500 m")
lines!(ax_F1000, 1:time_after, F1000[1,:], color = :red, label = "1000 m")
lines!(ax_F1000, 1:time_after, F2000[1,:], color = :blue, label = "2000 m")
axislegend(ax_F1000, position = :rt)

ax_F140 = Axis(fig_flux[2, 1]; xlabel = "time (days)", ylabel = "F(z)/F₁₄₀", yaxisposition = :left)
ylims!(ax_F140, 0, 1)
lines!(ax_F140, 1:time_after, F200[1,:]./F140[1,:], color = :black, label = "200 m/140 m")
lines!(ax_F140, 1:time_after, F200_martin[1,:], color = :black, linestyle = :dash, label = "200 m/140 m (Martin)")
lines!(ax_F140, 1:time_after, F500[1,:]./F140[1,:], color = :green, label = "500 m/140 m")
lines!(ax_F140, 1:time_after, F500_martin[1,:], color = :green, linestyle = :dash, label = "500 m/140 m (Martin)")
lines!(ax_F140, 1:time_after, F1000[1,:]./F140[1,:], color = :red, label = "1000 m/140 m")
lines!(ax_F140, 1:time_after, F1000_martin[1,:], color = :red, linestyle = :dash, label = "1000 m/140 m (Martin)")
lines!(ax_F140, 1:time_after, F2000[1,:]./F140[1,:], color = :blue, label = "2000 m/140 m")
lines!(ax_F140, 1:time_after, F2000_martin[1,:], color = :blue, linestyle = :dash, label = "2000 m/140 m (Martin)")
axislegend(ax_F140, position = :rt)

display(fig_flux)
=#
