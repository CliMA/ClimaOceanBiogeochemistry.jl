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
using Oceananigans.Advection: FluxFormAdvection

const Ny = 500 
const Nz = 100
const Ly = 15000kilometers   # m
const Lz = 2000           # m

arch = CPU()
# We use a two-dimensional grid, with a `Flat` `y`-direction:
grid = RectilinearGrid(arch,
                       size = (Ny, Nz),
                       y = (0, Ly),
                       z = (-Lz, 0),
                       topology=(Flat, Bounded, Bounded))

# Set streamfunction
deltaN = 200kilometers   # North downwelling width
deltaS = 500kilometers  # South upwelling width
deltaZ = -1000     # Vertical asymmetry
Ψᵢ(y, z)  = - 1.2 * ((1 - exp(-y / deltaS)) * (1 - exp(-(Ly - y) / deltaN)) * 
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

# (I have to specify u to allow CheckPointer)
u = XFaceField(grid) 
fill_halo_regions!(u, arch)

############################# Model setup ############################# 

kz(y,z,t) = 1e-4 + 5e-3 * (tanh((z+150)/20)+1) + 1e-2 * exp(-(z+2000)/50)
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
maxNCP(y,z) = (1e-5 + sinpi(y/Ly) * 5e-5)/day
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
                                                                        incident_PAR = incident_PAR,
                                                                        particulate_organic_phosphorus_sedremin_timescale = 0.5 / day),
                                                                        #iron_scavenging_rate = 0),
                                    velocities = PrescribedVelocityFields(; u, v, w),
                                    tracers = (:DIC, :ALK, :PO₄, :NO₃, :DOP, :POP, :Fe),
                                    # forcing=(; Fe=F_forcing),
                                    tracer_advection = FluxFormAdvection(
                                                                    Centered(),
                                                                    Centered(),
                                                                    Centered()),
                                    coriolis = nothing,
                                    buoyancy = nothing,
                                    closure = (tracer_vertical_closure, tracer_horizontal_closure))

set!(model, DIC=2.1, ALK=2.35, NO₃=2e-2, PO₄=1.2e-3, DOP=0, POP=0, Fe = 6e-7) # mol PO₄ m⁻³

simulation = Simulation(model; Δt = 1days, stop_time=300days) 

# Define a callback to zero out Fe tendency
function modify_tendency!(model)
        model.timestepper.Gⁿ.Fe .= 0
        return nothing
end                                        
simulation.callbacks[:modify_Fe] = Callback(modify_tendency!, 
                                            callsite = TendencyCallsite())

# Print the progress 
progress(sim) = @printf("Iteration: %d, time: %s\n", 
                        iteration(sim), prettytime(sim)) 
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
                        schedule = TimeInterval(10days), 
                        filename = "AMOC_test",
                        overwrite_existing = true)

# simulation.output_writers[:checkpointer] = Checkpointer(model2,
#             schedule = TimeInterval(180days),
#             prefix = "AMOC26_checkpoint",
#             overwrite_existing = false)
        
run!(simulation) # , pickup = false)

#################################### Visualize ####################################

filepath = simulation.output_writers[:simple_output].filepath
# filepath = "./AMOC64.jld2"

v_timeseries = FieldTimeSeries(filepath, "v")
w_timeseries = FieldTimeSeries(filepath, "w")

times = w_timeseries.times
xw, yw, zw = nodes(w_timeseries)

PO4_timeseries = FieldTimeSeries(filepath, "PO₄")
# times = PO4_timeseries.times
# xw, yw, zw = nodes(PO4_timeseries)
POP_timeseries = FieldTimeSeries(filepath, "POP")
DOP_timeseries = FieldTimeSeries(filepath, "DOP")
Fe_timeseries = FieldTimeSeries(filepath, "Fe")
NO3_timeseries = FieldTimeSeries(filepath, "NO₃")

n = Observable(1)
# title = @lift @sprintf("t = Year %d x 50", times[$n] / (3652.5*5)days) 
title = @lift @sprintf("t = Day %d0", times[$n] / 10days) 

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
colors6 = cgrad(:RdYlBu_6, 6, categorical=true, rev=true) #Reverse(:RdYlBu_6)
colors8 = cgrad(:RdYlBu_8, 8, categorical=true, rev=true) #Reverse(:RdYlBu_8)

ax_s = Axis(fig[2, 1]; xlabel = "y (km)", ylabel = "z (m)",title = "Stream function (m² s⁻¹)", aspect = 1)
hm_s = heatmap!(ax_s, yw./1e3, zw, Ψ; colormap = :viridis) 
Colorbar(fig[2, 2], hm_s; flipaxis = false)
contour!(ax_s, yw./1e3, grid.zᵃᵃᶜ[0:100], Ψ, levels = 10, color = :black)

ax_v = Axis(fig[2, 3]; xlabel = "y (km)", ylabel = "z (m)",title = "v (t) (cm s⁻¹)", aspect = 1)
hm_v = heatmap!(ax_v, yw./1e3, zw, vₙ; colorrange = (-2.5,2.5), colormap = :balance) 
Colorbar(fig[2, 4], hm_v; flipaxis = false)

ax_w = Axis(fig[2, 5]; xlabel = "y (km)", ylabel = "z (m)", title = "w (t) (cm s⁻¹)", aspect = 1)
hm_w = heatmap!(ax_w, yw./1e3, zw, wₙ; colorrange = (-8e-4, 8e-4), colormap = :balance) 
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
hm_Fe = heatmap!(ax_Fe, yw/1e3, zw, Feₙ; colorrange = (0,1.5),colormap = colors6) 
Colorbar(fig[3, 6], hm_Fe; flipaxis = false)

ax_POP = Axis(fig[4, 3]; xlabel = "y (km)", ylabel = "z (m)", title = "[POP] (μM)", aspect = 1)
hm_POP = heatmap!(ax_POP, yw/1e3, zw, POPₙ; colorrange = (0,0.01),colormap = :rainbow1) 
Colorbar(fig[4, 4], hm_POP; flipaxis = false)

ax_avg_PO4 = Axis(fig[4, 1:2]; xlabel = "[PO₄] (μM)", ylabel = "z (m)", title = "Average [PO₄] (μM)", yaxisposition = :right)
xlims!(ax_avg_PO4, 0, 4)
PO4_prof = lines!(ax_avg_PO4, avg_PO4ₙ[][1, :], zw[1:100])

ax_avg_POP = Axis(fig[4, 5:6]; xlabel = "[POP] (μM)", ylabel = "z (m)", title = "Average [POP] (μM)",yaxisposition = :right)
xlims!(ax_avg_POP, 0, 0.01)
POP_prof = lines!(ax_avg_POP, avg_POPₙ[][1, :], zw[1:100])

fig[1, 1:6] = Label(fig, title, tellwidth=false)

# And, finally, we record a movie.
frames = 1:length(times)
record(fig, "AMOC_test.mp4", frames, framerate=25) do i
    n[] = i
    PO4_prof[1] = avg_PO4ₙ[][1, :]
    POP_prof[1] = avg_POPₙ[][1, :]
end
nothing #hide

################################ Video comparing two models ################################
# filepath1 = "./AMOC64.jld2"
# filepath2 = "./AMOC63.jld2"

# PO4_timeseries1 = FieldTimeSeries(filepath1, "PO₄")
# PO4_timeseries2 = FieldTimeSeries(filepath2, "PO₄")
# times = PO4_timeseries1.times
# xw, yw, zw = nodes(PO4_timeseries1)

# n = Observable(1)
# title = @lift @sprintf("t = Year %d x 50", times[$n] / (3652.5*5)days) 

# PO4ₙ1 = @lift 1e3*interior(PO4_timeseries1[$n], 1, :, :)
# avg_PO4ₙ1 = @lift mean(1e3*interior(PO4_timeseries1[$n], 1, :, :), dims=1) 
# PO4ₙ2 = @lift 1e3*interior(PO4_timeseries2[$n], 1, :, :)
# avg_PO4ₙ2 = @lift mean(1e3*interior(PO4_timeseries2[$n], 1, :, :), dims=1) 

# fig = Figure(size = (800, 700))
# colors6 = cgrad(:RdYlBu_6, 6, categorical=true, rev=true) #Reverse(:RdYlBu_6)
# colors8 = cgrad(:RdYlBu_8, 8, categorical=true, rev=true) #Reverse(:RdYlBu_8)

# ax_PO4 = Axis(fig[2, 1]; xlabel = "y (km)", ylabel = "z (m)", title = "[PO₄] (μM)", aspect = 1)
# hm_PO4 = heatmap!(ax_PO4, yw/1e3, zw, PO4ₙ1; colorrange = (0,3),colormap = colors6, interpolate = false) 
# Colorbar(fig[2, 2], hm_PO4; flipaxis = false)

# ax_PO42 = Axis(fig[3, 1]; xlabel = "y (km)", ylabel = "z (m)", title = "[PO₄] (μM)", aspect = 1)
# hm_PO42 = heatmap!(ax_PO42, yw/1e3, zw, PO4ₙ2; colorrange = (0,3),colormap = colors6, interpolate = false) 
# Colorbar(fig[3, 2], hm_PO42; flipaxis = false)

# ax_avg_PO4 = Axis(fig[2, 3]; xlabel = "[PO₄] (μM)", ylabel = "z (m)", title = "Average [PO₄] (μM)")
# xlims!(ax_avg_PO4, 0, 4)
# PO4_prof = lines!(ax_avg_PO4, avg_PO4ₙ1[][1, :], zw[1:100])

# ax_avg_PO42 = Axis(fig[3, 3]; xlabel = "[PO₄] (μM)", ylabel = "z (m)", title = "Average [PO₄] (μM)")
# xlims!(ax_avg_PO42, 0, 4)
# PO4_prof2 = lines!(ax_avg_PO42, avg_PO4ₙ2[][1, :], zw[1:100])

# fig[1, 1:3] = Label(fig, title, tellwidth=false)

# frames = 1:length(times)
# record(fig, "AMOC63vs64.mp4", frames, framerate=25) do i
#     n[] = i
#     PO4_prof[1] = avg_PO4ₙ1[][1, :]
#     PO4_prof2[1] = avg_PO4ₙ2[][1, :]
# end
# nothing #hide

# POP_timeseries1 = FieldTimeSeries(filepath1, "POP")
# POP_timeseries2 = FieldTimeSeries(filepath2, "POP")
# DOP_timeseries1 = FieldTimeSeries(filepath1, "DOP")
# DOP_timeseries2 = FieldTimeSeries(filepath2, "DOP")
# PO4_final1 = 1e3*interior(PO4_timeseries1[end], 1, :, :)
# PO4_final2 = 1e3*interior(PO4_timeseries2[end], 1, :, :)
# POP_final1 = 1e3*interior(POP_timeseries1[end], 1, :, :)
# POP_final2 = 1e3*interior(POP_timeseries2[end], 1, :, :)
# DOP_final1 = 1e3*interior(DOP_timeseries1[end], 1, :, :)
# DOP_final2 = 1e3*interior(DOP_timeseries2[end], 1, :, :)

# totalP1 = sum((PO4_final1 + POP_final1 + DOP_final1)*3e4*20)/1e12
# totalP2 = sum((PO4_final2 + POP_final2 + DOP_final2)*3e4*20)/1e12

# PO4_init = 1e3*interior(PO4_timeseries2[1], 1, :, :)
# DOP_init = 1e3*interior(DOP_timeseries2[1], 1, :, :)
# POP_init = 1e3*interior(POP_timeseries2[1], 1, :, :)
# totalP_init = sum((PO4_init + POP_init + DOP_init)*3e4*20)/1e12

####################################################################################
#=
PO4_final = 1e3*interior(PO4_timeseries[end], 1, :, :)
NO3_final = 1e3*interior(NO3_timeseries[end],1, :, :)
POP_final = 1e3*interior(POP_timeseries[end], 1, :, :)
DOP_final = 1e3*interior(DOP_timeseries[end], 1, :, :)
Fe_final = 1e6*interior(Fe_timeseries[end], 1, :, :)

z_matrix = repeat(zw[1:100], 1, 500)
# μᵖ= 3e-5 # /day
μ(y,z)= 1e-5+ sinpi(y/Ly) * 5e-5
μ_flat = zeros(length(yw), length(zw[1:100]))  
yw_broadcast = reshape(yw, length(yw), 1) 
zw_broadcast = reshape(zw[1:100], 1, length(zw[1:100])) 
μ_flat .= μ.(yw_broadcast, zw_broadcast)  

kᴵ=10
kᴾ=1e-7*1024.5 * 1e3 # unit conversion 1e3
kᴺ=1.6e-6*1024.5 * 1e3 # unit conversion 1e3
kᶠ=1e-10*1024.5 * 1e6 # unit conversion 1e6

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
min_lim_vals = min.(p_lim, n_lim, f_lim)  
NCP_final= max.(0, μ_flat .* light_lim .* min_lim_vals)
# NCP_result = μᵖ .* light_lim .* min_lim_vals 
NCP_result .= NCP_final .*1e3 # mmol P m⁻³ d⁻¹

z₀ = log(0.01)*25 
rₛₑ = 0.5 #/day
Remin_final = -0.84 .* 10 ./(z_matrix' .+ z₀) .* POP_final # mmol P m⁻³ d⁻¹
# Remin_final = ifelse(z < -1990, rₛₑ .* POP_final, 0.1 .* POP_final)
Remin_final[:,1] = rₛₑ .* POP_final[:,1] # sedimentary remin

# tracer concentration and limitation terms
fig_can = Figure(size = (1000, 800))

# Light limitation
# ax_l_lim = Axis(fig_can[1, 1]; xlabel = "y (km)", ylabel = "z (m)",title = "I/(I+kᵢ)", aspect = 1)
# hm_l_lim = heatmap!(ax_l_lim, yw/1e3, zw, light_lim; colorrange = (0, 1), colormap = :tempo) 
# Colorbar(fig_can[1, 2], hm_l_lim; flipaxis = false)
# Zoom into surface
ax_l_lim2 = Axis(fig_can[1, 1]; xlabel = "y (km)", ylabel = "z (m)",title = "Surface I/(I+kᵢ)", aspect = 1)
hm_l_lim2 = heatmap!(ax_l_lim2, yw/1e3, zw[85:100], light_lim[:,85:100]; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_can[1, 2], hm_l_lim2; flipaxis = false)
# NO3 limitation
ax_n_lim = Axis(fig_can[2, 1]; xlabel = "y (km)", ylabel = "z (m)",title = "NO₃/(NO₃+kₙₒ₃)", aspect = 1)
hm_n_lim = heatmap!(ax_n_lim, yw/1e3, zw, n_lim; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_can[2, 2], hm_n_lim; flipaxis = false)
# PO4 limitation
ax_p_lim = Axis(fig_can[2, 3]; xlabel = "y (km)", ylabel = "z (m)",title = "PO₄/(PO₄+kₚₒ₄)", aspect = 1)
hm_p_lim = heatmap!(ax_p_lim, yw/1e3, zw, p_lim; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_can[2, 4], hm_p_lim; flipaxis = false)
# Fe limitation
ax_f_lim = Axis(fig_can[1, 3]; xlabel = "y (km)", ylabel = "z (m)",title = "Fe/(Fe+kᵢᵣₒₙ)", aspect = 1)
hm_f_lim = heatmap!(ax_f_lim, yw/1e3, zw, f_lim; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_can[1, 4], hm_f_lim; flipaxis = false)
# Fe limitation
# ax_f_lim2 = Axis(fig_can[3, 3]; xlabel = "y (km)", ylabel = "z (m)",title = "Surface Fe/(Fe+kᵢᵣₒₙ)", aspect = 1)
# hm_f_lim2 = heatmap!(ax_f_lim2, yw/1e3, zw[85:100], f_lim[:,85:100]; colorrange = (0, 1), colormap = :tempo) 
# Colorbar(fig_can[3, 4], hm_f_lim2; flipaxis = false)

display(fig_can)

# PO4 concentrations VS WOA23
fig_can2 = Figure(size = (600, 1000))
data = PO4_final
ax_PO4 = Axis(fig_can2[1, 1]; xlabel = "y (km)", ylabel = "z (m)", title = "[PO₄] (μM)", aspect = 1)
hm_PO4 = heatmap!(ax_PO4, yw/1e3, zw, data; colorrange = (0, 3),
                colormap = colors6, interpolate = false)                 
Colorbar(fig_can2[1, 2], hm_PO4; flipaxis = false)
levels = range(0, 3, length=7) 
contour!(ax_PO4, yw./1e3, zw[1:100], data, levels = levels, color = :black)

data2 = NO3_final
ax_NO3 = Axis(fig_can2[2, 1]; xlabel = "y (km)", ylabel = "z (m)", title = "[NO₃] (μM)", aspect = 1)
hm_NO3 = heatmap!(ax_NO3, yw/1e3, zw, data2; colorrange = (0, 40),
                colormap = colors8, interpolate = false)                 
Colorbar(fig_can2[2, 2], hm_NO3; flipaxis = false)
levels2 = range(0, 40, length=9) 
contour!(ax_NO3, yw./1e3, zw[1:100], data2, levels = levels2, color = :black)
display(fig_can2)

#=
Fe_final = 1e6*interior(Fe_timeseries[end], 1, :, :)
F = max.(0, Fe_final)
kᶠ=1e-10*1024.5 * 1e6 
f_lim = F ./ (F .+ kᶠ)

fig_can2 = Figure(size = (800, 400))

ax_Fe = Axis(fig_can2[1, 1]; xlabel = "y (km)", ylabel = "z (m)", title = "[Fe] (nM)", aspect = 1)
hm_Fe = heatmap!(ax_Fe, yw/1e3, zw, Fe_final; colorrange = (0,1),colormap = Reverse(:RdYlBu_6)) 
Colorbar(fig_can2[1, 2], hm_Fe; flipaxis = false)
contour!(ax_Fe, yw./1e3, zw[1:100], Fe_final, levels = 5, color = :black)

ax_f_lim = Axis(fig_can2[1, 3]; xlabel = "y (km)", ylabel = "z (m)",title = "Fe/(Fe+kᵢᵣₒₙ)", aspect = 1)
hm_f_lim = heatmap!(ax_f_lim, yw/1e3, zw, f_lim; colorrange = (0, 1), colormap = :tempo) 
Colorbar(fig_can2[1, 4], hm_f_lim; flipaxis = false)
contour!(ax_f_lim, yw/1e3, zw[1:100], f_lim, levels = 5, color = :black)
display(fig_can2)
=#
################################## Plot prod/remin rates ##################################

fig_rate = Figure(size = (800, 1200))

ax_NCP = Axis(fig_rate[1, 1]; xlabel = "y (km)", ylabel = "z (m)",title = "Top NCP (mmol m⁻³ d⁻¹)", aspect = 1)
hm_NCP = heatmap!(ax_NCP, yw/1e3, zw[85:100], NCP_result[:,85:100]; colorrange = (0,4e-2),colormap = :viridis) 
Colorbar(fig_rate[1, 2], hm_NCP; flipaxis = false)

ax_remin = Axis(fig_rate[1, 3]; xlabel = "y (km)", ylabel = "z (m)",title = "POP remin (mmol m⁻³ d⁻¹)", aspect = 1)
hm_remin = heatmap!(ax_remin, yw/1e3, zw, Remin_final; colorrange = (0,1.2e-3),colormap = :viridis) 
Colorbar(fig_rate[1, 4], hm_remin; flipaxis = false)

avg_NCP = mean(NCP_result; dims = 1)
ax_avg_NCP = Axis(fig_rate[2, 1:2]; xlabel = "Rate (mmol m⁻³ d⁻¹)", ylabel = "z (m)", title = "Ave. NCP")
xlims!(ax_avg_NCP, -5e-4, 2.5e-2)
lines!(ax_avg_NCP, vec(avg_NCP),zw[1:100], linewidth = 2)

avg_R = mean(Remin_final; dims = 1)
ax_avg_R = Axis(fig_rate[2, 3:4]; xlabel = "Rate (mmol m⁻³ d⁻¹)", ylabel = "z (m)", title = "Ave. POP Remin")
xlims!(ax_avg_R, -1e-6, 1.5e-3)
lines!(ax_avg_R, vec(avg_R),zw[1:100], linewidth = 2)

total_NCP = sum(NCP_result*20; dims = 2) # mmol m⁻² d⁻¹
ax_total_NCP = Axis(fig_rate[3, 1:2]; xlabel = "y (km)", ylabel = "Rate (mmol m⁻² d⁻¹)", title = "Integrated NCP")
# xlims!(ax_avg_NCP, -1e-4, 0.5)
lines!(ax_total_NCP, yw/1e3, vec(total_NCP),linewidth = 2)

total_R = sum((Remin_final + DOP_final./30).*20; dims = 2)
ax_total_R = Axis(fig_rate[3, 3:4]; xlabel = "y (km)", ylabel = "Rate (mmol m⁻² d⁻¹)", title = "Integrated (POP + DOP) Remin")
# xlims!(ax_avg_R, -1e-6, 2e-3)
lines!(ax_total_R, yw/1e3, vec(total_R),linewidth = 2)

zi = zw[1:100]
# compare to Martin curve
Remin_flux = vec(mean(10 .* POP_final; dims = 1))
ref_grid_index = 91
Martin_flux = Remin_flux[ref_grid_index]*((zi[ref_grid_index]+z₀)./(zi.+z₀)).^0.84
Martin_flux[ref_grid_index:end] .= NaN

ax_flux = Axis(fig_rate[4, 3:4]; xlabel = "POP flux (mmol m⁻² d⁻¹)", ylabel = "z (m)", title = "Flux comparison")
xlims!(ax_flux, 0, 0.2)
lines!(ax_flux, Remin_flux[2:100], zi[2:100], linewidth = 3, label = "Model")
lines!(ax_flux, Martin_flux, zi, linewidth = 2, label = "Martin curve")
axislegend(ax_flux, position = :rb)

# Calculate e-ratio 

eratio_120 = 10 .* POP_final[:,95] ./ total_NCP # mmol m⁻² d⁻¹

ax_eR = Axis(fig_rate[4, 1:2]; xlabel = "y (km)", ylabel = "e-ratio", title = "e-ratio (120 m)")
ylims!(ax_eR, 0, 0.2)
# Remove boundaries: 0 PAR = 0 NCP
lines!(ax_eR, yw[3:498]/1e3, vec(eratio_120[3:498]), linewidth = 2)

display(fig_rate)

=#
