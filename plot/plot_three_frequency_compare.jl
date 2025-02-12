using GLMakie
using Printf
using Statistics

using Oceananigans
using Oceananigans.Units
using Oceananigans.Fields: ZeroField, CenterField, FunctionField
using Oceananigans.BoundaryConditions: fill_halo_regions!

filepath1 = "./AMOC115_1perY.jld2"
filepath2 = "./AMOC115_4perY.jld2"
filepath3 = "./AMOC115_12perY.jld2"

# Load POP data
POP_timeseries1 = FieldTimeSeries(filepath1, "POP")
times = POP_timeseries1.times
xw, yw, zw = nodes(POP_timeseries1)
POP_timeseries2 = FieldTimeSeries(filepath2, "POP")
POP_timeseries3 = FieldTimeSeries(filepath3, "POP")

############################## vertical diffusivity ##############################
kz1(y,z,t) = 1e-4 + 5e-3 * (tanh((z+(100+50*sinpi(2*(t/day-(365.25*2000 + 1))/(365))))/20)+1) + 1e-2 * exp(-(z+4000)/50)
kz2(y,z,t) = 1e-4 + 5e-3 * (tanh((z+(100+50*sinpi(2*(t/day-(365.25*2000 + 1))/(365/4))))/20)+1) + 1e-2 * exp(-(z+4000)/50)
kz3(y,z,t) = 1e-4 + 5e-3 * (tanh((z+(100+50*sinpi(2*(t/day-(365.25*2000 + 1))/(365/12))))/20)+1) + 1e-2 * exp(-(z+4000)/50)

t= collect(1days:1days:366days)
kz1_prof = Matrix{Float64}(undef, length(zw), length(t))
kz2_prof = Matrix{Float64}(undef, length(zw), length(t))
kz3_prof = Matrix{Float64}(undef, length(zw), length(t))
for (i, t_val) in enumerate(t)
    kz1_prof[:, i] = [kz1(0, z_val, t_val) for z_val in zw]
    kz2_prof[:, i] = [kz2(0, z_val, t_val) for z_val in zw]
    kz3_prof[:, i] = [kz3(0, z_val, t_val) for z_val in zw]
end

n = Observable(1)
title = @lift @sprintf("t = Day %d", ((times[$n] - 365.25*2000days) / 1days)-1) 

####################### POP flux = (10m/d)*[POP] #######################
# 1D profile
avg_POP_flux_init1 = mean(1e4*interior(POP_timeseries1[1], 1, :, :) , dims=1) 
avg_POP_flux1ₙ = @lift mean(1e4*interior(POP_timeseries1[$n], 1, :, :), dims=1) 

avg_POP_flux_init2 = mean(1e4*interior(POP_timeseries2[1], 1, :, :) , dims=1) 
avg_POP_flux2ₙ = @lift mean(1e4*interior(POP_timeseries2[$n], 1, :, :), dims=1) 

avg_POP_flux_init3 = mean(1e4*interior(POP_timeseries3[1], 1, :, :) , dims=1) 
avg_POP_flux3ₙ = @lift mean(1e4*interior(POP_timeseries3[$n], 1, :, :), dims=1) 

# 2D heatmap
POP_flux_init1 = 1e4*interior(POP_timeseries1[1], 1, :, :) 
POP_flux1ₙ = @lift 1e4*interior(POP_timeseries1[$n], 1, :, :)

POP_flux_init2 = 1e4*interior(POP_timeseries2[1], 1, :, :) 
POP_flux2ₙ = @lift 1e4*interior(POP_timeseries2[$n], 1, :, :)

POP_flux_init3 = 1e4*interior(POP_timeseries3[1], 1, :, :) 
POP_flux3ₙ = @lift 1e4*interior(POP_timeseries3[$n], 1, :, :)

################## Martin curve ##################
z₀ = log(0.01)*25 
Martin_factor = ((zw[190] .+ z₀) ./ (zw .+ z₀)).^0.84

# 1D martin
avg_Martin_flux1ₙ = @lift($avg_POP_flux1ₙ[190] .* Martin_factor)
avg_Martin_flux2ₙ = @lift($avg_POP_flux2ₙ[190] .* Martin_factor)
avg_Martin_flux3ₙ = @lift($avg_POP_flux3ₙ[190] .* Martin_factor)

# 2D heatmap
Martin_flux1ₙ = Observable(similar(POP_flux1ₙ[]))  
matrix1 = Observable(similar(POP_flux1ₙ[]))
Martin_flux1ₙ = @lift begin
    for i in 1:200
        $matrix1[:, i] .= $POP_flux1ₙ[:, 190] .* $Martin_factor[i]
    end
    return matrix1
end
diff_Martin1 = @lift ($POP_flux1ₙ .- $Martin_flux1ₙ[])

Martin_flux2ₙ = Observable(similar(POP_flux2ₙ[]))  
matrix2 = Observable(similar(POP_flux2ₙ[]))
Martin_flux2ₙ = @lift begin
    for i in 1:200
        $matrix2[:, i] .= $POP_flux2ₙ[:, 190] .* $Martin_factor[i]
    end
    return matrix2
end
diff_Martin2 = @lift ($POP_flux2ₙ .- $Martin_flux2ₙ[])

Martin_flux3ₙ = Observable(similar(POP_flux3ₙ[])) 
matrix3 = Observable(similar(POP_flux3ₙ[]))
Martin_flux3ₙ = @lift begin
    for i in 1:200
        $matrix3[:, i] .= $POP_flux3ₙ[:, 190] .* $Martin_factor[i]
    end
    return matrix3
end
diff_Martin3 = @lift ($POP_flux3ₙ .- $Martin_flux3ₙ[])

############################# Plot #############################
fig_compare = Figure(size=(1000, 1000))

##################
# 1 cycle / year
##################
ax_kz1 = Axis(fig_compare[2, 1:2], xlabel = "kz (m² s⁻¹)", ylabel = "z (m)", title = "kz")
kz_prof1 = lines!(ax_kz1, kz1_prof[:, 1], zw)
ylims!(ax_kz1, -1000, 0)

# POP flux 
ax_avg_POP_flux1 = Axis(fig_compare[3, 1:2]; xlabel = "POP flux (mmol m⁻² d⁻¹)", ylabel = "z (m)", title = "Average POP flux")
lines!(ax_avg_POP_flux1, vec(avg_POP_flux_init1), zw, color = :grey, linewidth = 1, label = "Initial")
POP_flux_prof1 = lines!(ax_avg_POP_flux1, avg_POP_flux1ₙ[][1,:], zw,  linewidth = 2, label = "Dynamic")
Martin_flux_prof1 = lines!(ax_avg_POP_flux1, avg_Martin_flux1ₙ[], zw, color = :red,  linewidth = 2, linestyle = :dash, label = "Martin")
ylims!(ax_avg_POP_flux1, -1000, 0)
xlims!(ax_avg_POP_flux1, 0, 0.08)
axislegend(ax_avg_POP_flux1, position = :rb)

# diff POP flux - Martin
ax_diffMartin1 = Axis(fig_compare[4, 1]; xlabel = "y (km)", ylabel = "z (m)", title = "Δ(POP flux): True - Martin", aspect = 1)
hm_diffMartin1 = heatmap!(ax_diffMartin1, yw/1e3, zw, diff_Martin1; colorrange = (-7e-3,7e-3),colormap = :balance) 
Colorbar(fig_compare[4, 2], hm_diffMartin1; flipaxis = false)
ylims!(ax_diffMartin1, -1000, -200)

##################
# 4 cycles / year
##################
ax_kz2 = Axis(fig_compare[2, 3:4], xlabel = "kz (m² s⁻¹)", ylabel = "z (m)", title = "kz")
kz_prof2 = lines!(ax_kz2, kz2_prof[:, 1], zw)
ylims!(ax_kz2, -1000, 0)

# POP flux 
ax_avg_POP_flux2 = Axis(fig_compare[3, 3:4]; xlabel = "POP flux (mmol m⁻² d⁻¹)", ylabel = "z (m)", title = "Average POP flux")
lines!(ax_avg_POP_flux2, vec(avg_POP_flux_init2), zw, color = :grey, linewidth = 1, label = "Initial")
POP_flux_prof2 = lines!(ax_avg_POP_flux2, avg_POP_flux2ₙ[][1,:], zw,  linewidth = 2, label = "Dynamic")
Martin_flux_prof2 = lines!(ax_avg_POP_flux2, avg_Martin_flux2ₙ[], zw, color = :red,  linewidth = 2, linestyle = :dash, label = "Martin")
ylims!(ax_avg_POP_flux2, -1000, 0)
xlims!(ax_avg_POP_flux2, 0, 0.08)
axislegend(ax_avg_POP_flux2, position = :rb)

# diff POP flux - Martin
ax_diffMartin2 = Axis(fig_compare[4, 3]; xlabel = "y (km)", ylabel = "z (m)", title = "Δ(POP flux): True - Martin", aspect = 1)
hm_diffMartin2 = heatmap!(ax_diffMartin2, yw/1e3, zw, diff_Martin2; colorrange = (-7e-3,7e-3),colormap = :balance) 
Colorbar(fig_compare[4, 4], hm_diffMartin2; flipaxis = false)
ylims!(ax_diffMartin2, -1000, -200)

##################
# 12 cycles / year
##################
ax_kz3 = Axis(fig_compare[2, 5:6], xlabel = "kz (m² s⁻¹)", ylabel = "z (m)", title = "kz")
kz_prof3 = lines!(ax_kz3, kz3_prof[:, 1], zw)
ylims!(ax_kz3, -1000, 0)

# POP flux 
ax_avg_POP_flux3 = Axis(fig_compare[3, 5:6]; xlabel = "POP flux (mmol m⁻² d⁻¹)", ylabel = "z (m)", title = "Average POP flux")
lines!(ax_avg_POP_flux3, vec(avg_POP_flux_init3), zw, color = :grey, linewidth = 1, label = "Initial")
POP_flux_prof3 = lines!(ax_avg_POP_flux3, avg_POP_flux3ₙ[][1,:], zw,  linewidth = 2, label = "Dynamic")
Martin_flux_prof3 = lines!(ax_avg_POP_flux3, avg_Martin_flux3ₙ[], zw, color = :red,  linewidth = 2, linestyle = :dash, label = "Martin")
ylims!(ax_avg_POP_flux3, -1000, 0)
xlims!(ax_avg_POP_flux3, 0, 0.08)
axislegend(ax_avg_POP_flux3, position = :rb)

# diff POP flux - Martin
ax_diffMartin3 = Axis(fig_compare[4, 5]; xlabel = "y (km)", ylabel = "z (m)", title = "Δ(POP flux): True - Martin", aspect = 1)
hm_diffMartin3 = heatmap!(ax_diffMartin3, yw/1e3, zw, diff_Martin3; colorrange = (-7e-3,7e-3),colormap = :balance) 
Colorbar(fig_compare[4, 6], hm_diffMartin3; flipaxis = false)
ylims!(ax_diffMartin3, -1000, -200)

fig_compare[1, 1:6] = Label(fig_compare, title, tellwidth=false)

# And, finally, we record a movie.
frames = 1:length(times)
record(fig_compare, "AMOC115_3flux_compare.mp4", frames, framerate=50) do i
    n[] = i
    kz_prof1[1]= kz1_prof[:,i]
    kz_prof2[1]= kz2_prof[:,i]
    kz_prof3[1]= kz3_prof[:,i]
    POP_flux_prof1[1] = avg_POP_flux1ₙ[][1,:]
    POP_flux_prof2[1] = avg_POP_flux2ₙ[][1,:]
    POP_flux_prof3[1] = avg_POP_flux3ₙ[][1,:]
    Martin_flux_prof1[1] = avg_Martin_flux1ₙ[]
    Martin_flux_prof2[1] = avg_Martin_flux2ₙ[]
    Martin_flux_prof3[1] = avg_Martin_flux3ₙ[]
end
nothing #hide