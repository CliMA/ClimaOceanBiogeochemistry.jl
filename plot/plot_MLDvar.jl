using GLMakie
using Printf
using Statistics

ds = range(0, 365, length=365*2)  # Time in days
MLD_1 = 100 .+ 50 .* sinpi.(2 .* (ds ./ 365)) 
# MLD_2 = MLD_1 .+ 40 .* sinpi.(2 .* (ds ./ (365/12))) 
# MLD_3 = MLD_1 .+ 40 .* sinpi.(2 .* (ds ./ (365/52))) 

PAR_1 = 700 .* (1 .+ sinpi.((ds ./ 365)))
PAR_2 = 700 .* sinpi(15/15000) .* (1 .+ sinpi.((ds ./ 365)))

fig = Figure(size=(800, 400))

ax1 = Axis(fig[1, 1], xlabel="Day of the Year", ylabel="MLD (m)", title="Seasonal (yearly)")
lines!(ax1, ds, MLD_1,color=:grey)
ylims!(ax1,200,0)

ax = Axis(fig[1, 2], xlabel="Day of the Year", ylabel="Equator PAR (W m⁻²)", title="Seasonal (yearly)")
lines!(ax, ds, PAR_1,color=:grey)
ylims!(ax,700,1400)

ax_low = Axis(fig[1, 2]; xlabel="Day of the Year", ylabel="Polar PAR (W m⁻²)",  yaxisposition=:right, 
                yticks=0:1:5, yticklabelcolor=:red3, ylabelcolor=:red3)
lines!(ax_low, ds, PAR_2,color=:red3)
ylims!(ax_low, 2,5)
linkxaxes!(ax, ax_low)

fig

# ax = Axis(fig[2, 1], xlabel="Day of the Year", ylabel="MLD (m)", title="Seasonal (yearly) + mesoscale eddy (monthly)")
# lines!(ax, ds, MLD_2,color=:grey)
# ylims!(ax,200,0)

# ax = Axis(fig[3, 1], xlabel="Day of the Year", ylabel="MLD (m)", title="Seasonal (yearly) + submesoscale eddy (weekly)")
# lines!(ax, ds, MLD_3,color=:grey)
# ylims!(ax,200,0)
# fig
