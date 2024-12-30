using GLMakie
# using CUDA
using Printf
using Statistics

using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients
using Oceananigans
using Oceananigans.Models.HydrostaticFreeSurfaceModels:
                    HydrostaticFreeSurfaceModel,
                    PrescribedVelocityFields
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization
using Oceananigans.Units
using Oceananigans.Fields: ZeroField
using Oceananigans.BoundaryConditions: fill_halo_regions!

Ny = 500 
Nz = 100
Ly = 15000kilometers   # m
Lz = 2000           # m

arch = CPU()
# We use a two-dimensional grid, with a `Flat` `y`-direction:
grid = RectilinearGrid(arch,
                       size = (Ny, Nz),
                       y = (0, Ly),
                       z = (-Lz, 0),
                       topology=(Flat, Bounded, Bounded))

# Set streamfunction
deltaN = 500kilometers   # North downwelling width
deltaS = 1000kilometers  # South upwelling width
deltaZ = 1000     # Vertical asymmetry
Ψᵢ(y, z)  = -8 * ((1 - exp(-y / deltaS)) * (1 - exp(-(Ly - y) / deltaN)) * 
            sinpi(z / Lz) * exp(z/deltaZ))

# Ψᵢ(y, z) = 10*(y/Ly) * sinpi(-z/Lz) * exp((2/Lz)*z) * (1-exp((y-Ly)/(0.2*Ly)))
# Ψᵢ(x, z) = 0.1 * sinpi(x/Lx) * sinpi(-z/Lz) * exp((2/Lz)*z) * exp((2/Lx)*x)
Ψ = Field{Face, Center, Face}(grid)
set!(Ψ, Ψᵢ)
fill_halo_regions!(Ψ, arch)

# Set velocity field from streamfunction
v = YFaceField(grid)
w = ZFaceField(grid)
v .= - ∂z(Ψ)
w .= + ∂y(Ψ)
fill_halo_regions!(v, arch)
fill_halo_regions!(w, arch)

# (I have to specify v to allow CheckPointer)
u = XFaceField(grid) 
fill_halo_regions!(u, arch)

############################# Model setup ############################# 

kz(y,z,t) = 5e-4 + 5e-3 * (tanh((z+150)/20)+1) + 1e-2 * exp(-(z+2000)/50)
vertical_closure = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), κ=kz)
horizontal_closure = HorizontalScalarDiffusivity(κ=1e3)

# Set PAR as a function of latitude
incident_PAR = CenterField(grid) 
surface_PAR(y,z) = 700 * sinpi(y/Ly) 
set!(incident_PAR, surface_PAR)   
fill_halo_regions!(incident_PAR, arch)

# Model
model = HydrostaticFreeSurfaceModel(grid = grid,
                                    biogeochemistry = CarbonAlkalinityNutrients(; grid,
                                                                        maximum_net_community_production_rate  = 1e-5/day,
                                                                        incident_PAR = incident_PAR,
                                                                        particulate_organic_phosphorus_sedremin_timescale = 0.5 / day,
                                                                        iron_scavenging_rate = 0),
                                    velocities = PrescribedVelocityFields(; u, v, w),
                                    tracers = (:DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe),
                                    coriolis = nothing,
                                    buoyancy = nothing,
                                    closure = (vertical_closure, horizontal_closure))

set!(model, DIC=2.1, ALK=2.35, NO₃=2e-2, PO₄=1.2e-3, DOP=0, POP=0, Fe = 5e-7) # mol PO₄ m⁻³

simulation = Simulation(model; Δt = 1days, stop_time=1000days) 

# We add a `TimeStepWizard` callback to adapt the simulation's time-step,
# wizard = TimeStepWizard(cfl=0.2, max_change=1.1, max_Δt=60minutes)
# simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(50))

# Print the progress 
progress(sim) = @printf("Iteration: %d, time: %s\n", 
                        iteration(sim), prettytime(sim)) 
add_callback!(simulation, progress, IterationInterval(2000))

outputs = (v = model.velocities.v,
            w = model.velocities.w,
            PO₄= model.tracers.PO₄,
            DOP = model.tracers.DOP,
            POP = model.tracers.POP,
            NO₃ = model.tracers.NO₃,
            Fe = model.tracers.Fe)

simulation.output_writers[:simple_output] =
        JLD2OutputWriter(model, outputs, 
                        schedule = TimeInterval(1days), 
                        filename = "AMOC19",
                        overwrite_existing = true)
        
# simulation.output_writers[:checkpointer] = Checkpointer(model,
#                         schedule = TimeInterval(365.25*500days),
#                         prefix = "AMOC14_checkpoint",
#                         overwrite_existing = true)
        
run!(simulation) #, pickup = false)
        
#################################### Visualize ####################################

filepath = simulation.output_writers[:simple_output].filepath
# filepath = "./AMOC17.jld2"

v_timeseries = FieldTimeSeries(filepath, "v")
w_timeseries = FieldTimeSeries(filepath, "w")

times = w_timeseries.times
xw, yw, zw = nodes(w_timeseries)

PO4_timeseries = FieldTimeSeries(filepath, "PO₄")
POP_timeseries = FieldTimeSeries(filepath, "POP")
DOP_timeseries = FieldTimeSeries(filepath, "DOP")
Fe_timeseries = FieldTimeSeries(filepath, "Fe")
NO3_timeseries = FieldTimeSeries(filepath, "NO₃")

n = Observable(1)
# title = @lift @sprintf("t = Year %d0", times[$n] / 3652.5days) 
title = @lift @sprintf("t = Day %d", times[$n] / 1days) 

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

fig = Figure(size = (1500, 1500))

ax_s = Axis(fig[2, 1]; xlabel = "y (km)", ylabel = "z (m)",title = "Stream function (m² s⁻¹)", aspect = 1)
hm_s = heatmap!(ax_s, yw./1e3, zw, Ψ; colormap = :viridis) 
Colorbar(fig[2, 2], hm_s; flipaxis = false)
contour!(ax_s, yw./1e3, grid.zᵃᵃᶜ[0:100], Ψ, levels = 10, color = :black)

ax_v = Axis(fig[2, 3]; xlabel = "y (km)", ylabel = "z (m)",title = "v (t) (cm s⁻¹)", aspect = 1)
hm_v = heatmap!(ax_v, yw./1e3, zw, vₙ; colorrange = (-2,2), colormap = :balance) 
Colorbar(fig[2, 4], hm_v; flipaxis = false)

ax_w = Axis(fig[2, 5]; xlabel = "y (km)", ylabel = "z (m)", title = "w (t) (cm s⁻¹)", aspect = 1)
hm_w = heatmap!(ax_w, yw./1e3, zw, wₙ; colorrange = (-5e-4, 5e-4), colormap = :balance) 
Colorbar(fig[2, 6], hm_w; flipaxis = false)

ax_PO4 = Axis(fig[3, 1]; xlabel = "y (km)", ylabel = "z (m)", title = "[PO₄] (μM)", aspect = 1)
hm_PO4 = heatmap!(ax_PO4, yw/1e3, zw, PO4ₙ; colorrange = (0,2),colormap = :jet1) 
Colorbar(fig[3, 2], hm_PO4; flipaxis = false)

ax_NO3 = Axis(fig[3, 3]; xlabel = "y (km)", ylabel = "z (m)", title = "[NO₃] (μM)", aspect = 1)
hm_NO3 = heatmap!(ax_NO3, yw/1e3, zw, NO3ₙ; colorrange = (0,30),colormap = :jet1) 
Colorbar(fig[3, 4], hm_NO3; flipaxis = false)

ax_Fe = Axis(fig[3, 5]; xlabel = "y (km)", ylabel = "z (m)", title = "[Fe] (nM)", aspect = 1)
hm_Fe = heatmap!(ax_Fe, yw/1e3, zw, Feₙ; colorrange = (0,0.8),colormap = :jet1) 
Colorbar(fig[3, 6], hm_Fe; flipaxis = false)

ax_POP = Axis(fig[4, 3]; xlabel = "y (km)", ylabel = "z (m)", title = "[POP] (μM)", aspect = 1)
hm_POP = heatmap!(ax_POP, yw/1e3, zw, POPₙ; colorrange = (0,0.03),colormap = :jet1) 
Colorbar(fig[4, 4], hm_POP; flipaxis = false)

ax_avg_PO4 = Axis(fig[4, 1:2]; xlabel = "[PO₄] (μM)", ylabel = "z (m)", title = "Average [PO₄] (μM)", yaxisposition = :right)
xlims!(ax_avg_PO4, 0, 2)
PO4_prof = lines!(ax_avg_PO4, avg_PO4ₙ[][1, :], zw[1:100])

ax_avg_POP = Axis(fig[4, 5:6]; xlabel = "[POP] (μM)", ylabel = "z (m)", title = "Average [POP] (μM)",yaxisposition = :right)
xlims!(ax_avg_POP, 0, 0.02)
POP_prof = lines!(ax_avg_POP, avg_POPₙ[][1, :], zw[1:100])

fig[1, 1:6] = Label(fig, title, tellwidth=false)

# And, finally, we record a movie.
frames = 1:length(times)
record(fig, "AMOC19.mp4", frames, framerate=50) do i
    n[] = i
    PO4_prof[1] = avg_PO4ₙ[][1, :]
    POP_prof[1] = avg_POPₙ[][1, :]
end
nothing #hide

####################################################################################
#=
PO4_final = 1e3*interior(PO4_timeseries[end], 1, :, :)
POP_final = 1e3*interior(POP_timeseries[end], 1, :, :)
DOP_final = 1e3*interior(DOP_timeseries[end], 1, :, :)
Fe_final = 1e6*interior(Fe_timeseries[end], 1, :, :)
NO3_final = 1e3*interior(NO3_timeseries[end],1, :, :)

z_matrix = repeat(zw[1:100], 1, 500)
μᵖ= 1e-5 # /day
kᴵ=10
kᴾ=1e-7*1024.5 * 1e3 # unit conversion 1e3
kᴺ=1.6e-6*1024.5 * 1e3 # unit conversion 1e3
kᶠ=1e-10*1024.5 * 1e6 # unit conversion 1e6
I = max.(0, 700 .* exp.(z_matrix' ./ 25))
P = max.(0, PO4_final)
N = max.(0, NO3_final)
F = max.(0, Fe_final)

light_lim = I ./ (I .+ kᴵ)
p_lim = P ./ (P .+ kᴾ)
n_lim = N ./ (N .+ kᴺ)
f_lim = F ./ (F .+ kᶠ)

NCP_final = max.(0, μᵖ .* light_lim .* min.(p_lim, n_lim, f_lim))
   
z₀ = log(0.01)*25 
rₛₑ = 0.5 #/day
Remin_final = -0.84 .* 10 ./(z_matrix' .+ z₀) .* POP_final
# Remin_final = ifelse(z < -1990, rₛₑ .* POP_final, 0.1 .* POP_final)
Remin_final[:,1:2] = rₛₑ .* POP_final[:,1:2] # sedimentary remin

# tracer concentration and limitation terms
fig_can = Figure(size = (800, 800))

# Light limitation
ax_l_lim = Axis(fig_can[1, 1]; xlabel = "y (km)", ylabel = "z (m)",title = "I/(I+kᵢ)", aspect = 1)
hm_l_lim = heatmap!(ax_l_lim, yw/1e3, zw, light_lim; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_can[1, 2], hm_l_lim; flipaxis = false)

ax_n_lim = Axis(fig_can[1, 3]; xlabel = "y (km)", ylabel = "z (m)",title = "NO₃/(NO₃+kₙₒ₃)", aspect = 1)
hm_n_lim = heatmap!(ax_n_lim, yw/1e3, zw, n_lim; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_can[1, 4], hm_n_lim; flipaxis = false)

ax_p_lim = Axis(fig_can[2, 1]; xlabel = "y (km)", ylabel = "z (m)",title = "PO₄/(PO₄+kₚₒ₄)", aspect = 1)
hm_p_lim = heatmap!(ax_p_lim, yw/1e3, zw, p_lim; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_can[2, 2], hm_p_lim; flipaxis = false)

ax_f_lim = Axis(fig_can[2, 3]; xlabel = "y (km)", ylabel = "z (m)",title = "Fe/(Fe+kᵢᵣₒₙ)", aspect = 1)
hm_f_lim = heatmap!(ax_f_lim, yw/1e3, zw, f_lim; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_can[2, 4], hm_f_lim; flipaxis = false)

display(fig_can)
=#

# PO4 concentrations
# ax_PO4 = Axis(fig_can[1, 3]; xlabel = "y (km)", ylabel = "z (m)", title = "[PO₄] (μM)", aspect = 1)
# hm_PO4 = heatmap!(ax_PO4, yw/1e3, zw, PO4_final; colorrange = (0, 3.1),colormap = :jet1) 
# Colorbar(fig_can[1, 4], hm_PO4; flipaxis = false)
# # PO4 limitation
# ax_p_lim = Axis(fig_can[1, 5]; xlabel = "y (km)", ylabel = "z (m)",title = "PO₄/(PO₄+kₚₒ₄)", aspect = 1)
# hm_p_lim = heatmap!(ax_p_lim, yw/1e3, zw, p_lim; colorrange = (0, 1), colormap = :tempo) 
# Colorbar(fig_can[1, 6], hm_p_lim; flipaxis = false)
# POP concentrations
# ax_POP = Axis(fig_can[2, 1]; xlabel = "y (km)", ylabel = "z (m)", title = "[POP] (μM)", aspect = 1)
# hm_POP = heatmap!(ax_POP, yw/1e3, zw, POP_final; colorrange = (0,0.012),colormap = :jet1) 
# Colorbar(fig_can[2, 2], hm_POP; flipaxis = false)
# # NO3 concentrations
# ax_NO3 = Axis(fig_can[2, 3]; xlabel = "y (km)", ylabel = "z (m)", title = "[NO₃] (μM)", aspect = 1)
# hm_NO3 = heatmap!(ax_NO3, yw/1e3, zw, NO3_final; colorrange = (0,50), colormap = :jet1) 
# Colorbar(fig_can[2, 4], hm_NO3; flipaxis = false)
# NO3 limitation
# ax_n_lim = Axis(fig_can[2, 5]; xlabel = "y (km)", ylabel = "z (m)",title = "NO₃/(NO₃+kₙₒ₃)", aspect = 1)
# hm_n_lim = heatmap!(ax_n_lim, yw/1e3, zw, n_lim; colorrange = (0, 1), colormap = :tempo) 
# Colorbar(fig_can[2, 6], hm_n_lim; flipaxis = false)
# DOP concentrations
# ax_DOP = Axis(fig_can[3, 1]; xlabel = "y (km)", ylabel = "z (m)", title = "[DOP] (μM)", aspect = 1)
# hm_DOP = heatmap!(ax_DOP, yw/1e3, zw, DOP_final; colorrange = (0,0.06),colormap = :jet1) 
# Colorbar(fig_can[3, 2], hm_DOP; flipaxis = false)
# # Fe concentrations
# ax_Fe = Axis(fig_can[3, 3]; xlabel = "y (km)", ylabel = "z (m)", title = "[Fe] (nM)", aspect = 1)
# hm_Fe = heatmap!(ax_Fe, yw/1e3, zw, Fe_final; colorrange = (0,1.5), colormap = :jet1) 
# Colorbar(fig_can[3, 4], hm_Fe; flipaxis = false)
# Fe limitation
# ax_f_lim = Axis(fig_can[3, 5]; xlabel = "y (km)", ylabel = "z (m)",title = "Fe/(Fe+kᵢᵣₒₙ)", aspect = 1)
# hm_f_lim = heatmap!(ax_f_lim, yw/1e3, zw, f_lim; colorrange = (0, 1), colormap = :tempo) 
# Colorbar(fig_can[3, 6], hm_f_lim; flipaxis = false)

# display(fig_can)

################################## Plot prod/remin rates ##################################
#=
fig_rate = Figure(size = (1000, 1000))

ax_NCP = Axis(fig_rate[1, 1]; xlabel = "y (km)", ylabel = "z (m)",title = "NCP (mmol m⁻³ d⁻¹)", aspect = 1)
hm_NCP = heatmap!(ax_NCP, yw/1e3, zw, 1000*NCP_final; colorrange = (0,1e-2),colormap = :jet1) 
Colorbar(fig_rate[1, 2], hm_NCP; flipaxis = false)

fig_rate = Figure(size = (500, 500))
ax_NCP = Axis(fig_rate[1, 1]; xlabel = "y (km)", ylabel = "z (m)",title = "Top NCP (mmol m⁻³ d⁻¹)", aspect = 1)
hm_NCP = heatmap!(ax_NCP, yw/1e3, zw[85:100], 1000*NCP_final[:,85:100]; colorrange = (0,1e-2),colormap = :jet1) 
Colorbar(fig_rate[1, 2], hm_NCP; flipaxis = false)
display(fig_rate)

ax_remin = Axis(fig_rate[1, 3]; xlabel = "y (km)", ylabel = "z (m)",title = "Remin (mmol m⁻³ d⁻¹)", aspect = 1)
hm_remin = heatmap!(ax_remin, yw/1e3, zw, Remin_final; colorrange = (0,5e-4),colormap = :jet1) 
Colorbar(fig_rate[1, 4], hm_remin; flipaxis = false)

avg_NCP = mean(1000*NCP_final; dims = 1)
ax_avg_NCP = Axis(fig_rate[2, 1:2]; xlabel = "Rate (mmol m⁻³ d⁻¹)", ylabel = "z (m)", title = "Ave. NCP", yaxisposition = :right)
xlims!(ax_avg_NCP, -1e-4, 1e-2)
lines!(ax_avg_NCP, vec(avg_NCP),zw[1:100])

avg_R = mean(Remin_final; dims = 1)
ax_avg_R = Axis(fig_rate[2, 3:4]; xlabel = "Rate (mmol m⁻³ d⁻¹)", ylabel = "z (m)", title = "Ave. POP Remin", yaxisposition = :right)
xlims!(ax_avg_R, -1e-6, 5e-4)
lines!(ax_avg_R, vec(avg_R),zw[1:100])

zi = zw[1:100]
# compare to Martin curve
Remin_flux = vec(mean(10 .* POP_final; dims = 1))
Martin_flux = Remin_flux[92]*((zi[92]+z₀)./(zi.+z₀)).^0.84

ax_flux = Axis(fig_rate[3, 3:4]; xlabel = "POP flux (mmol m⁻² d⁻¹)", ylabel = "z (m)", title = "Flux comparison", yaxisposition = :right)
xlims!(ax_flux, 0, 0.08)
lines!(ax_flux, Remin_flux[3:100], zi[3:100], label = "Model")
lines!(ax_flux, Martin_flux, zi, label = "Martin curve")
axislegend(ax_flux, position = :rb)

display(fig_rate)
=#