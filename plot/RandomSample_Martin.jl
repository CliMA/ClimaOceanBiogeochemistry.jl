using GLMakie
using Printf
using Statistics

using Oceananigans
using Oceananigans.Units

#################################################################
##### Load data: Premin (rate) and POP concentration (flux) #####
#################################################################

filepath1 = "./1perY_11y.jld2"
filepath2 = "./4perY_11y.jld2"
filepath3 = "./12perY_11y.jld2"

POP_timeseries1 = FieldTimeSeries(filepath1, "POP")
times = POP_timeseries1.times
xw, yw, zw = nodes(POP_timeseries1)

POP_timeseries2 = FieldTimeSeries(filepath2, "POP")
POP_timeseries3 = FieldTimeSeries(filepath3, "POP")

#################################################################
######## Process data: flux (t,y,z) ########
#################################################################

POP_flux1 = zeros(365, 500, 200) 
POP_flux2 = zeros(365, 500, 200) 
POP_flux3 = zeros(365, 500, 200) 
for i in 1:365
    POP_flux1[i, :, :] = 1e4*interior(POP_timeseries1[i], 1, :, :)
    POP_flux2[i, :, :] = 1e4*interior(POP_timeseries2[i], 1, :, :)
    POP_flux3[i, :, :] = 1e4*interior(POP_timeseries3[i], 1, :, :)
end


z₀ = log(0.01)*25 
Martin_factor = ((zw[190] .+ z₀) ./ (zw .+ z₀)).^0.84


#################################################################
####### Random sampling: any time (t), any column (y) #######
#################################################################

# Number of random samples
num_samples = 20

# Generate random indices
t_ind1 = rand(1:365, num_samples)
y_ind1 = rand(1:500, num_samples)
t_ind2 = rand(1:365, num_samples)
y_ind2 = rand(1:500, num_samples) 
t_ind3 = rand(1:365, num_samples)
y_ind3 = rand(1:500, num_samples)

rand_flux1 = zeros(num_samples, 200)
rand_flux2 = zeros(num_samples, 200) 
rand_flux3 = zeros(num_samples, 200) 
for ind in 1:num_samples
    rand_flux1[ind,:] = POP_flux1[t_ind1[ind],y_ind1[ind],:]
    rand_flux2[ind,:] = POP_flux2[t_ind2[ind],y_ind2[ind],:]
    rand_flux3[ind,:] = POP_flux3[t_ind3[ind],y_ind3[ind],:]
end

avg_flux1 = mean(rand_flux1; dims = 1)
avg_flux2 = mean(rand_flux2; dims = 1)
avg_flux3 = mean(rand_flux3; dims = 1)

Martin_flux1 = avg_flux1[190].* Martin_factor
Martin_flux2 = avg_flux2[190].* Martin_factor
Martin_flux3 = avg_flux3[190].* Martin_factor

#################################################################
######################### Plot fluxes #########################
#################################################################

fig = Figure(size=(1200, 400))

ax1 = Axis(fig[1, 1]; ylabel = "z (m)", xlabel = "flux (mmol m⁻² d⁻¹)", title = "POP flux")
for ind in 1:num_samples
    scatter!(ax1, rand_flux1[ind,:], zw, color=:grey, markersize=5)
end
lines!(ax1, vec(avg_flux1), zw, color=:black, linewidth = 3, label = "Mean")
lines!(ax1, Martin_flux1, zw, color=:red3, linewidth = 2, label = "Martin")
ylims!(ax1, -1000, 0)
axislegend(ax1, position = :rb)

ax2 = Axis(fig[1, 2]; ylabel = "z (m)", xlabel = "flux (mmol m⁻² d⁻¹)", title = "POP flux")
for ind in 1:num_samples
    scatter!(ax2, rand_flux2[ind,:], zw, color=:grey, markersize=5)
end
lines!(ax2, vec(avg_flux2), zw, color=:black, linewidth = 3, label = "Mean")
lines!(ax2, Martin_flux2, zw, color=:red3, linewidth = 2, label = "Martin")
ylims!(ax2, -1000, 0)
axislegend(ax2, position = :rb)

ax3 = Axis(fig[1, 3]; ylabel = "z (m)", xlabel = "flux (mmol m⁻² d⁻¹)", title = "POP flux")
for ind in 1:num_samples
    scatter!(ax3, rand_flux3[ind,:], zw, color=:grey, markersize=5)
end
lines!(ax3, vec(avg_flux3), zw, color=:black, linewidth = 3, label = "Mean")
lines!(ax3, Martin_flux3, zw, color=:red3, linewidth = 2, label = "Martin")
ylims!(ax3, -1000, 0)
axislegend(ax3, position = :rb)

display(fig)