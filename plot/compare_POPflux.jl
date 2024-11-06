using JLD2
using GLMakie
using Oceananigans
using Statistics

#= Mouw (2016) POP flux compiled: 
obs_POPflux = JLD2.jldopen("plot/POPflux_Mouw2016.jld2", "r") do file
    file["obs_POPflux"]
end

# Note that the original POP flux has a unit of mg/m2/day
# Convert mg to mmol: no clear value. 
# Use 3553 g/mol here as an example
obs_POPflux[:,2] = obs_POPflux[:,2] ./ 3553
=#

obs_POPflux = JLD2.jldopen("plot/POCflux_Cael2018.jld2", "r") do file
    file["obs_POCflux"]
end

# Note that the original POC flux has a unit of mg C/m2/d
# Convert POC to POP (mmol/m2/day): follow Redfield C:P value. 
obs_POPflux[:,2] = obs_POCflux[:,2] ./ 106 ./12.01

# Load two models
POP_constR_timeseries = FieldTimeSeries("CAN_2D_constR_0d1.jld2", "POP")
POP_constR_final = 1e3*interior(POP_constR_timeseries[end], :, 1, :)

# POP_Rz_timeseries = FieldTimeSeries("./CAN_2D_Rz7.jld2", "POP")
# POP_Rz_final = 1e3*interior(POP_Rz_timeseries[end], :, 1, :)

times = POP_constR_timeseries.times
xw, yw, zw = nodes(POP_constR_timeseries)

# Rz_flux = vec(mean(10 .* POP_Rz_final; dims = 1))
constR_flux = vec(mean(10 .* POP_constR_final; dims = 1))

z₀ = log(0.01)*25 
Martin_flux = constR_flux[190]*(zw./z₀).^-0.84
# Martin_flux = Rz_flux[183]*((zw[183]+z₀)./(zw.+z₀)).^0.84

fig = Figure(size = (500, 500))
ax_flux = Axis(fig[1, 1]; xlabel = "POP flux (mmol m⁻² d⁻¹)", ylabel = "z (m)", title = "POP flux comparisons", yaxisposition = :right)

xlims!(ax_flux, 0, 0.55)
ylims!(ax_flux, -2000, 0)

# Observational data
sc_flux = scatter!(ax_flux, obs_POPflux[:,2],-obs_POPflux[:,1], marker = :circle, color = :black, transparency = true, alpha = 0.5,label = "Observations (Cael et al., 2018)")

# lines!(ax_flux, Rz_flux, zw, label = "r = bxwₛ/(z+z₀)")
lines!(ax_flux, constR_flux, zw, label = "r = 0.018/day")
lines!(ax_flux, Martin_flux, zw, label = "Martin curve")

axislegend(ax_flux, position = :rb)

display(fig)