using GLMakie
using Oceananigans
using Statistics

z₀ = -115.0 #log(0.01)*25 
z₁ = -1010.0
z₂ = -210.0
Δz₁ = z₁-z₀
Δz₂ = z₂-z₀
martin_z = ((z₁+z₀)/(z₂+z₀))^-0.84

w₁ = -10 #m/d
w₂ = -50

t_range = collect(1:0.01:365)
T= collect(1:1:365) # period

range_Teff1 = zeros(1, 365) 
range_Teff2 = zeros(1, 365)

amplitude = 0.25

for i in 1:365
    ratio1 = martin_z .* (1.0 .+ amplitude .* sinpi.(2/T[i].*(t_range.- Δz₁/w₁))) ./ (1.0 .+ amplitude .* sinpi.(2/T[i].*(t_range.-Δz₂/w₁)))
    range_Teff1[i] = maximum(ratio1)-minimum(ratio1) #std(ratio)
    ratio2 = martin_z .* (1.0 .+ amplitude .* sinpi.(2/T[i].*(t_range.- Δz₁/w₂))) ./ (1.0 .+ amplitude .* sinpi.(2/T[i].*(t_range.-Δz₂/w₂)))
    range_Teff2[i] = maximum(ratio2)-minimum(ratio2) #std(ratio)
end

ratio_T3 = martin_z .* (1.0 .+ amplitude .* sinpi.(2/T[3].*(t_range.- Δz₁/w₁))) ./ (1.0 .+ amplitude .* sinpi.(2/T[3].*(t_range.-Δz₂/w₁)))
ratio_T10 = martin_z .* (1.0 .+ amplitude .* sinpi.(2/T[10].*(t_range.- Δz₁/w₁))) ./ (1.0 .+ amplitude .* sinpi.(2/T[10].*(t_range.-Δz₂/w₁)))
# ratio_T30 = martin_z .* (1.0 .+ amplitude .* sinpi.(2/T[30].*(t_range.- Δz₁/w₁))) ./ (1.0 .+ amplitude .* sinpi.(2/T[30].*(t_range.-Δz₂/w₁)))
# ratio_T160 = martin_z .* (1.0 .+ amplitude .* sinpi.(2/T[160].*(t_range.- Δz₁/w₁))) ./ (1.0 .+ amplitude .* sinpi.(2/T[160].*(t_range.-Δz₂/w₁)))
# ratio_T365 = martin_z .* (1.0 .+ amplitude .* sinpi.(2/T[365].*(t_range.- Δz₁/w₁))) ./ (1.0 .+ amplitude .* sinpi.(2/T[365].*(t_range.-Δz₂/w₁)))


fig = Figure()

ax_tT1 = Axis(fig[1, 1]; ylabel = "F₁₀₀₀/F₂₀₀", xlabel = "t (day)", title = "F₁₀₀₀/F₂₀₀ vs time (w = 10 m/d)")
lines!(ax_tT1, t_range, ratio_T3, linewidth = 2, label = "T = 3 d")
lines!(ax_tT1, t_range, ratio_T10, linewidth = 2, label = "T = 10 d")
# lines!(ax_tT1, t_range, ratio_T30, linewidth = 2, label = "T = 30 d")
# lines!(ax_tT1, t_range, ratio_T160, linewidth = 2, label = "T = 160 d")
# lines!(ax_tT1, t_range, ratio_T365, linewidth = 2, label = "T = 365 d")
axislegend(ax_tT1, position = :rt)

fig = Figure()
ax_tT2 = Axis(fig[2, 1]; ylabel = "ΔF₁₀₀₀/F₂₀₀ (max-min)", xlabel = "Period T (day)", title = "ΔF₁₀₀₀/F₂₀₀ vs Period T",
            yticklabelcolor=:brown3, ylabelcolor=:brown3)
lines!(ax_tT2, T, vec(range_Teff1), color=:brown3, linewidth = 2.5, label = "w = 10 m/d")
lines!(ax_tT2, T, vec(range_Teff2), color=:black, linewidth = 2.5, label = "w = 50 m/d")
xlims!(ax_tT2, 0,10)
axislegend(ax_tT2, position = :lt)

display(fig)