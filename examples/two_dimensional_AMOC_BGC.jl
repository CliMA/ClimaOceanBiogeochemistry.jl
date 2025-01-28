using GLMakie
# using CUDA
using Printf
using Statistics

using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients
using Oceananigans
using Oceananigans.Units
using Oceananigans.Fields: ZeroField, CenterField
using Oceananigans.BoundaryConditions: fill_halo_regions!

using Oceananigans.Models.HydrostaticFreeSurfaceModels:
                    HydrostaticFreeSurfaceModel,
                    PrescribedVelocityFields
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization
using Oceananigans: TendencyCallsite
# using Oceananigans.Advection: FluxFormAdvection

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

# Set streamfunction
deltaN = 1000kilometers   # North downwelling width
deltaS = 2000kilometers  # South upwelling width
deltaZ = 1000     # Vertical asymmetry
Ψᵢ(y, z)  = - 12 * ((1 - exp(-y / deltaS)) * (1 - exp(-(Ly - y) / deltaN)) * 
            sinpi(z / Lz) * exp(z/deltaZ))

# Ψᵢ(y, z) = 10*(y/Ly) * sinpi(-z/Lz) * exp((2/Lz)*z) * (1-exp((y-Ly)/(0.2*Ly)))
# Ψᵢ(x, z) = 0.1 * sinpi(x/Lx) * sinpi(-z/Lz) * exp((2/Lz)*z) * exp((2/Lx)*x)
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

############################# Model setup ############################# 

kz(y,z,t) = 1e-4 + 5e-3 * (tanh((z+100)/20)+1) + 1e-2 * exp(-(z+4000)/50)
tracer_vertical_closure = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), 
                                                κ=(DIC=kz,ALK=kz,PO₄=kz,NO₃=kz,DOP=kz,POP=kz,Fe=0))
tracer_horizontal_closure = HorizontalScalarDiffusivity(κ=(DIC=1e3,ALK=1e3,PO₄=1e3,NO₃=1e3,DOP=1e3,POP=1e3,Fe=0))

# Set Fe as a frozen field
# @inline Feᵢ(y,z) = (0.02+(-z/2000)) /1e6

# Fe dust flux
# Fe_input(y, z, take!) = (z > -30 ? 1e-15 * sinpi(y/Ly) : 0.0)
# F_forcing = Forcing(Fe_input)

# Set max NCP as a function of latitude (y)
maximum_net_community_production_rate = CenterField(grid)
maxNCP(y,z) = (8e-6 + sinpi(y/Ly) * 8e-5)/day
set!(maximum_net_community_production_rate, maxNCP)
fill_halo_regions!(maximum_net_community_production_rate, arch)

# Set PAR as a function of latitude
incident_PAR = CenterField(grid) 
surface_PAR(y,z) = 700 * sinpi(y/Ly) 
set!(incident_PAR, surface_PAR)   
fill_halo_regions!(incident_PAR, arch)

# Model
model = HydrostaticFreeSurfaceModel(grid = grid,
                                    biogeochemistry = CarbonAlkalinityNutrients(; grid,
                                                                        maximum_net_community_production_rate  = maximum_net_community_production_rate,
                                                                        incident_PAR = incident_PAR),
                                                                        #particulate_organic_phosphorus_sedremin_timescale = 0.1 / day),
                                                                        #iron_scavenging_rate = 0),
                                    velocities = PrescribedVelocityFields(; u, v, w),
                                    tracers = (:DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe),
                                    # forcing=(; Fe=F_forcing),
                                    coriolis = nothing,
                                    buoyancy = nothing,
                                    closure = (tracer_vertical_closure, tracer_horizontal_closure))

set!(model, DIC=2.1, ALK=2.35, NO₃=2.4e-2, PO₄=1.7e-3, DOP=0, POP=0, Fe = 6e-7) # mol PO₄ m⁻³

simulation = Simulation(model; Δt = 1day, stop_time=365.25days) 

# Define a callback to zero out Fe tendency
function modify_tendency!(model)
        model.timestepper.Gⁿ.Fe .= 0
        return nothing
end                                        
simulation.callbacks[:modify_Fe] = Callback(modify_tendency!, 
                                            callsite = TendencyCallsite())

# Print the progress 
progress(sim) = @printf("Iteration: %d, time: %s, total(P): %.2e \n",
            iteration(sim), prettytime(sim),
            sum(model.tracers.PO₄) + sum(model.tracers.POP) + sum(model.tracers.DOP))
add_callback!(simulation, progress, IterationInterval(100))

outputs = (v = model.velocities.v,
            w = model.velocities.w,
            PO₄= model.tracers.PO₄,
            DOP = model.tracers.DOP,
            POP = model.tracers.POP,
            NO₃ = model.tracers.NO₃,
            Fe = model.tracers.Fe)

simulation.output_writers[:simple_output] =
        JLD2OutputWriter(model, outputs, 
                        schedule = TimeInterval(365.25days), 
                        filename = "AMOC_test4",
                        overwrite_existing = true)

# simulation.output_writers[:checkpointer] = Checkpointer(model2,
#             schedule = TimeInterval(180days),
#             prefix = "AMOC26_checkpoint",
#             overwrite_existing = false)
        
run!(simulation) # , pickup = false)

#################################### Visualize ####################################

filepath = simulation.output_writers[:simple_output].filepath
# filepath = "./AMOC91.jld2"

v_timeseries = FieldTimeSeries(filepath, "v")
w_timeseries = FieldTimeSeries(filepath, "w")

times = w_timeseries.times
xv, yv, zv = nodes(v_timeseries)
xw, yw, zw = nodes(w_timeseries)

PO4_timeseries = FieldTimeSeries(filepath, "PO₄")
xt, yt, zt = nodes(PO4_timeseries)
POP_timeseries = FieldTimeSeries(filepath, "POP")
DOP_timeseries = FieldTimeSeries(filepath, "DOP")
Fe_timeseries = FieldTimeSeries(filepath, "Fe")
NO3_timeseries = FieldTimeSeries(filepath, "NO₃")

n = Observable(1)
title = @lift @sprintf("t = Year %d", times[$n] / 365.25days) 
# title = @lift @sprintf("t = Day %d0", times[$n] / 10days) 

# convert unit from m/s to cm/s:
vₙ = @lift 100*interior(v_timeseries[$n], 1, :, :)
wₙ = @lift 100*interior(w_timeseries[$n], 1, :, :)

# convert unit from mol/m³ to μM: 1e3*interior(...)
PO4ₙ = @lift 1e3*interior(PO4_timeseries[$n], 1, :, :)
avg_PO4ₙ = @lift mean(1e3*interior(PO4_timeseries[$n], 1, :, :), dims=1) 

POPₙ = @lift 1e3*interior(POP_timeseries[$n], 1, :, :)
avg_POPₙ = @lift mean(1e3*interior(POP_timeseries[$n], 1, :, :), dims=1) 

NO3ₙ = @lift 1e3*interior(NO3_timeseries[$n], 1, :, :)
Feₙ = @lift 1e6*interior(Fe_timeseries[$n], 1, :, :)

#################################################################################

colors5 = cgrad(:RdYlBu, 5, categorical=true, rev=true) 
colors6 = cgrad(:RdYlBu, 6, categorical=true, rev=true) 
colors8 = cgrad(:RdYlBu, 8, categorical=true, rev=true)


fig = Figure(size = (1500, 1500))

ax_s = Axis(fig[2, 1]; xlabel = "y (km)", ylabel = "z (m)",title = "Stream function (m² s⁻¹)", aspect = 1)
hm_s = heatmap!(ax_s, yw./1e3, zw, Ψ; colormap = :viridis) 
Colorbar(fig[2, 2], hm_s; flipaxis = false)
contour!(ax_s, yv./1e3, zw, Ψ, levels = 10, color = :black)

ax_v = Axis(fig[2, 3]; xlabel = "y (km)", ylabel = "z (m)",title = "v (t) (cm s⁻¹)", aspect = 1)
hm_v = heatmap!(ax_v, yw./1e3, zw, vₙ; colorrange = (-2.2,2.2), colormap = :balance) 
Colorbar(fig[2, 4], hm_v; flipaxis = false)

ax_w = Axis(fig[2, 5]; xlabel = "y (km)", ylabel = "z (m)", title = "w (t) (cm s⁻¹)", aspect = 1)
hm_w = heatmap!(ax_w, yw./1e3, zw, wₙ; colorrange = (-6e-4, 6e-4), colormap = :balance) 
Colorbar(fig[2, 6], hm_w; flipaxis = false)

ax_PO4 = Axis(fig[3, 1]; xlabel = "y (km)", ylabel = "z (m)", title = "[PO₄] (μM)", aspect = 1)
hm_PO4 = heatmap!(ax_PO4, yw/1e3, zw, PO4ₙ; colorrange = (0,3),colormap = colors6, interpolate = false) 
Colorbar(fig[3, 2], hm_PO4; flipaxis = false)
# Note: Adding Contours often cause crashes!!!
# contour!(ax_PO4, yw/1e3, zw[1:100], PO4ₙ, levels = 5, color = :black)

ax_NO3 = Axis(fig[3, 3]; xlabel = "y (km)", ylabel = "z (m)", title = "[NO₃] (μM)", aspect = 1)
hm_NO3 = heatmap!(ax_NO3, yw/1e3, zw, NO3ₙ; colorrange = (0,40),colormap = colors8, interpolate = false) 
Colorbar(fig[3, 4], hm_NO3; flipaxis = false)

ax_Fe = Axis(fig[3, 5]; xlabel = "y (km)", ylabel = "z (m)", title = "[Fe] (nM)", aspect = 1)
hm_Fe = heatmap!(ax_Fe, yw/1e3, zw, Feₙ; colorrange = (0,1),colormap = colors6) 
Colorbar(fig[3, 6], hm_Fe; flipaxis = false)

ax_POP = Axis(fig[4, 3]; xlabel = "y (km)", ylabel = "z (m)", title = "[POP] (μM)", aspect = 1)
hm_POP = heatmap!(ax_POP, yw/1e3, zw, POPₙ; colorrange = (0,0.006),colormap = :rainbow1) 
Colorbar(fig[4, 4], hm_POP; flipaxis = false)

ax_avg_PO4 = Axis(fig[4, 1:2]; xlabel = "[PO₄] (μM)", ylabel = "z (m)", title = "Average [PO₄] (μM)", yaxisposition = :right)
xlims!(ax_avg_PO4, 0, 3)
PO4_prof = lines!(ax_avg_PO4, avg_PO4ₙ[][1, :], zt)

ax_avg_POP = Axis(fig[4, 5:6]; xlabel = "[POP] (μM)", ylabel = "z (m)", title = "Average [POP] (μM)",yaxisposition = :right)
xlims!(ax_avg_POP, 0, 0.005)
POP_prof = lines!(ax_avg_POP, avg_POPₙ[][1, :], zt)

fig[1, 1:6] = Label(fig, title, tellwidth=false)

# And, finally, we record a movie.
frames = 1:length(times)
record(fig, "AMOC_test4.mp4", frames, framerate=10) do i
    n[] = i
    PO4_prof[1] = avg_PO4ₙ[][1, :]
    POP_prof[1] = avg_POPₙ[][1, :]
end
nothing #hide


#= Comfirm steady-state
num_timesteps = 41 # number of time slices saved
PO4_sum_vector = zeros(Float64, num_timesteps)
DOP_sum_vector = zeros(Float64, num_timesteps)
POP_sum_vector = zeros(Float64, num_timesteps)
NO3_sum_vector = zeros(Float64, num_timesteps)
PO4_single_vector = zeros(Float64, num_timesteps)
POP_single_vector = zeros(Float64, num_timesteps)
DOP_single_vector = zeros(Float64, num_timesteps)
NO3_single_vector = zeros(Float64, num_timesteps)

for index in 1:num_timesteps
    # mol P m-3
    PO4_sum_vector[index] = sum(interior(PO4_timeseries[index], 1, :, :))
    DOP_sum_vector[index] = sum(interior(DOP_timeseries[index], 1, :, :))
    POP_sum_vector[index] = sum(interior(POP_timeseries[index], 1, :, :))
    NO3_sum_vector[index] = sum(interior(NO3_timeseries[index], 1, :, :))
    PO4_single_vector[index] = PO4_timeseries[1, 100, 99, index]
    DOP_single_vector[index] = DOP_timeseries[1, 100, 99, index]
    POP_single_vector[index] = POP_timeseries[1, 100, 99, index]
    NO3_single_vector[index] = NO3_timeseries[1, 100, 99, index]
end

# Create the time vector `t`
t = 0:50:2000 

# Plot the data using GLMakie
fig_sum = Figure()
ax = Axis(fig_sum[1, 1], xlabel="Time (years)", ylabel="mol m⁻³", title="Sum of tracer vs Time")
lines!(ax, t, PO4_sum_vector./100, label="PO₄/100", linewidth = 3)
lines!(ax, t, DOP_sum_vector, label="DOP", linewidth = 3)
lines!(ax, t, POP_sum_vector, label="POP", linewidth = 3)
lines!(ax, t, NO3_sum_vector./1600, label="NO₃/1600", linewidth = 3)
axislegend(ax, position = (0.9,0.5))

ax_single = Axis(fig_sum[1, 2], xlabel="Time (years)", ylabel="mol m⁻³", title="A single point of tracer vs Time")
lines!(ax_single, t, PO4_single_vector./100, label="PO₄/100", linewidth = 3)
lines!(ax_single, t, DOP_single_vector./20, label="DOP/20", linewidth = 3)
lines!(ax_single, t, POP_single_vector, label="POP", linewidth = 3)
lines!(ax_single, t, NO3_single_vector./1600, label="NO₃/1600", linewidth = 3)
axislegend(ax_single, position = :rt)

display(fig_sum)
=#

################################# Load final data ################################

# POP_snd = 1e3*interior(POP_timeseries[end-1], 1, :, :)
#=
PO4_final = interior(PO4_timeseries[end], 1, :, :)
NO3_final = interior(NO3_timeseries[end], 1, :, :)
POP_final = interior(POP_timeseries[end], 1, :, :)
DOP_final = interior(DOP_timeseries[end], 1, :, :)
Fe_final = interior(Fe_timeseries[end], 1, :, :)

###################### Plot model nutrients vs WOA23 data #########################

using NCDatasets
dsP = Dataset("data/woa23_all_p00_01.nc")
dsN = Dataset("data/woa23_all_n00_01.nc")
lat = dsP["lat"][20:150]
# lon = dsP["lon"][:]
depth = dsP["depth"][:]
PO4_obs = dsP["p_an"][155,20:150,:] # Objectively analyzed mean fields for moles_of_phosphate_per_unit_mass_in_sea_water
NO3_obs = dsN["n_an"][155,20:150,:] # 25W

function interpolate_column(col)
    valid_indices = findall(!ismissing, col)
    
    # If column is all missing, return missing
    if isempty(valid_indices)
        return fill(missing, length(col))
    end
    
    valid_values = collect(Float64, col[valid_indices])
    
    # Create interpolation with boundary conditions
    itp = linear_interpolation(valid_indices, valid_values, 
                             extrapolation_bc=Flat())  # Use Flat() extrapolation
    
    # Fill all indices
    result = Vector{Union{Float64, Missing}}(undef, length(col))
    for i in 1:length(col)
        # Handle edge cases
        if i < minimum(valid_indices)
            result[i] = valid_values[1]  # Use first valid value
        elseif i > maximum(valid_indices)
            result[i] = valid_values[end]  # Use last valid value
        else
            result[i] = itp(i)
        end
    end
    
    return result
end
PO4_interpolated = hcat([interpolate_column(PO4_obs[:,j]) for j in 1:size(PO4_obs,2)]...)
NO3_interpolated = hcat([interpolate_column(NO3_obs[:,j]) for j in 1:size(NO3_obs,2)]...)

fig_can2 = Figure(size = (1000, 1000))

ax_PO4 = Axis(fig_can2[1, 1]; xlabel = "y (km)", ylabel = "z (m)", title = "[PO₄] (μM)", aspect = 1)
hm_PO4 = heatmap!(ax_PO4, yt/1e3, zt, PO4_final.*1e3; colorrange = (0, 3),
                colormap = colors6, interpolate = true)                 
Colorbar(fig_can2[1, 2], hm_PO4; flipaxis = false)
levels = range(0, 3, length=7) 
contour!(ax_PO4, yt./1e3, zt, PO4_final.*1e3, levels = levels, color = :black)

ax_PO4_obs = Axis(fig_can2[2, 1]; xlabel = "Latitude", ylabel = "z (m)", title = "[PO₄] (μM)", aspect = 1)
hm_PO4_obs = heatmap!(ax_PO4_obs, lat, -depth, PO4_interpolated; colorrange = (0, 3),
                colormap = colors6, interpolate = true)                 
Colorbar(fig_can2[2, 2], hm_PO4_obs; flipaxis = false)
ylims!(ax_PO4_obs, -4000, 0)
contour!(ax_PO4_obs, lat, -depth, PO4_interpolated, levels = levels, color = :black)

ax_NO3 = Axis(fig_can2[1, 3]; xlabel = "y (km)", ylabel = "z (m)", title = "[NO₃] (μM)", aspect = 1)
hm_NO3 = heatmap!(ax_NO3, yt/1e3, zt, NO3_final.*1e3; colorrange = (0, 40),
                colormap = colors8, interpolate = true)                 
Colorbar(fig_can2[1, 4], hm_NO3; flipaxis = false)
levels2 = range(0, 40, length=9) 
contour!(ax_NO3, yt./1e3, zt, NO3_final.*1e3, levels = levels2, color = :black)

ax_NO3_obs = Axis(fig_can2[2, 3]; xlabel = "Latitude", ylabel = "z (m)", title = "[NO₃] (μM)", aspect = 1)
hm_NO3_obs = heatmap!(ax_NO3_obs, lat, -depth, NO3_interpolated; colorrange = (0, 40),
                colormap = colors8, interpolate = true)                 
Colorbar(fig_can2[2, 4], hm_NO3_obs; flipaxis = false)
ylims!(ax_NO3_obs, -4000, 0)
contour!(ax_NO3_obs, lat, -depth, NO3_interpolated, levels = levels2, color = :black)
display(fig_can2)
=#

# ax_DOP = Axis(fig_can2[2, 1]; xlabel = "y (km)", ylabel = "z (m)", title = "[DOP] (μM)", aspect = 1)
# hm_DOP = heatmap!(ax_DOP, yt/1e3, zt, DOP_final.*1e3; colorrange = (0, 0.25),
#                 colormap = colors5, interpolate = false)                 
# Colorbar(fig_can2[2, 2], hm_DOP; flipaxis = false)

# ax_POP = Axis(fig_can2[2, 3]; xlabel = "y (km)", ylabel = "z (m)", title = "[POP] (μM)", aspect = 1)
# hm_POP = heatmap!(ax_POP, yt/1e3, zt, POP_final.*1e3; colorrange = (0, 0.006),
#                 colormap = colors6, interpolate = false)                 
# Colorbar(fig_can2[2, 4], hm_POP; flipaxis = false)

########################################################################
#=
z_matrix = repeat(zt, 1, 500)
kᴵ=30
kᴾ=5e-7 * 1024.5  # mol P m⁻³ 
kᴺ=7e-6 * 1024.5 
kᶠ=1e-10 * 1024.5  

incident_PAR = CenterField(grid) 
surface_PAR(y,z) = 700 * sinpi(y/Ly) 
set!(incident_PAR, surface_PAR)   
fill_halo_regions!(incident_PAR, arch)
I = incident_PAR[1,1:500,1:100] .* exp.(z_matrix' ./ 25)

# I = 700 .* exp.(z_matrix' ./ 25)
P = max.(0, PO4_final)
N = max.(0, NO3_final)
F = max.(0, Fe_final)

light_lim = I ./ (I .+ kᴵ)
p_lim = P ./ (P .+ kᴾ)
n_lim = N ./ (N .+ kᴺ)
f_lim = F ./ (F .+ kᶠ)
all_lim = light_lim .* min.(p_lim, n_lim, f_lim)  

############### Calculate Remin rates ###############
z₀ = log(0.01)*25 
rₛₑ = 0.5 # /day
wₛₙₖ = 10  # m/day
Δz = 20  # Vertical grid spacing in meters
POP_remin_final = -0.84 .* wₛₙₖ ./(z_matrix' .+ z₀) .* POP_final .*1e3 # mmol P m⁻³ d⁻¹
POP_remin_final[:,1] = rₛₑ .* POP_final[:,1] .*1e3 # sedimentary remin in the bottom layer

DOP_remin_final = DOP_final.* (2/365.25) .*1e3 # mmol P m⁻³ d⁻¹
total_remin_final = POP_remin_final .+ DOP_remin_final
=#

########################### Compare POC fluxes with observations ###########################
#=
using XLSX, JLD2, DataFrames

@load "data/POCflux_Cael2018.jld2" obs_POCflux
obs_depth_Cael = obs_POCflux[:,1]
obs_POC_flux_Cael = obs_POCflux[:,2]./ 12 .* 365.25 / 1e3 # mg C m⁻² d⁻¹ to mol C m⁻² y⁻¹
@load "data/POCflux_Cael2018Martin.jld2" obs_POCflux
obs_depth_Martin = obs_POCflux[:,1]
obs_POC_flux_Martin = obs_POCflux[:,2]./ 12 .* 365.25 / 1e3 # mg C m⁻² d⁻¹ to mol C m⁻² y⁻¹
 
# compare to Martin curve
POP_flux = vec(mean(wₛₙₖ .* POP_final .*1e3; dims = 1))
ref_grid_index = 91
Martin_flux = POP_flux[ref_grid_index]*((zt[ref_grid_index]+z₀)./(zt.+z₀)).^0.84
Martin_flux[ref_grid_index+3:end] .= NaN
model_POC_flux = POP_flux .* 117 .* 365.25  / 1e3# mmol P m⁻² d⁻¹ to mol C m⁻² y⁻¹ 

fig_comp = Figure(size = (800, 500))

ax_flux = Axis(fig_comp[1, 1]; xlabel = "POP flux (mmol P m⁻² d⁻¹)", ylabel = "z (m)", title = "POP flux")
xlims!(ax_flux, 0, 0.05)
lines!(ax_flux, POP_flux[2:100], zt[2:100], linewidth = 4, label = "Model")
lines!(ax_flux, Martin_flux, zt, linewidth = 2, color = :red3, label = "Martin curve")
axislegend(ax_flux, position = :rb)

ax_Cael= Axis(fig_comp[1, 2], xlabel = "POC flux (mol C m⁻² y⁻¹)", ylabel = "Depth (m)", title = "POC flux")
scatter!(ax_Cael, obs_POC_flux_Cael, -obs_depth_Cael, marker = :circle, color = :grey, label = "Cael et al. (2018)")
scatter!(ax_Cael, obs_POC_flux_Martin, obs_depth_Martin, marker = :circle, color = :red3, label = "Martin et al. (1987)")
lines!(ax_Cael, model_POC_flux[2:100], zt[2:100], color = :royalblue3, linewidth = 3, label = "Model")
xlims!(ax_Cael, 0, 3)
ylims!(ax_Cael, -2000, 0)
axislegend(ax_Cael, position = :rb)

display(fig_comp)
=#
