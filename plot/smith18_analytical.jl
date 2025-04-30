# Smith Jr, K. L., Ruhl, H. A., Huffard, C. L., Messié, M., & Kahru, M. (2018). 
# Episodic organic carbon fluxes from surface ocean to abyssal depths during long-term monitoring in NE Pacific. 
using GLMakie
using Printf
using Statistics
using Interpolations
using XLSX

xlsx_data = XLSX.readxlsx("Smith2018_PNAS.xlsx")
data2016 = xlsx_data["201617"]

day_year = data2016[2:end, 2]  # Column B (day)
EF = data2016[2:end, 4]  # 100 m
POC3400 = data2016[2:end, 3] 

#################################################################
###################### Analytical solution ######################
#################################################################

function compute_D(t::Vector, D0::Vector, z::Float64, z0::Float64, b::Float64, w::Float64)
    t_grid = reshape(t, 1, :)   # (1, Nt)
    t_shifted = clamp.(t_grid .- (z .- z0) ./ w, minimum(t), maximum(t))
    interp_D0 = linear_interpolation(t, D0)
    F_shifted = interp_D0.(t_shifted)
    Martin_term = ((z .+ z0) ./ (2*z0)).^(-b)
    return F_shifted .* Martin_term
end

function pointwise_w_fit(t, D0, D_pre, z, z0, b)
    w_opt = similar(t)
    interp_D0 = linear_interpolation(t, D0)

    for i in eachindex(t)
        t_i = t[i]
        D0_fun = interp_D0

        obj(w_i) = begin
            t_shift = clamp(t_i - (z - z0)/w_i, minimum(t), maximum(t))
            D_pred = D0_fun(t_shift) * ((z + z0) / (2 * z0))^(-b)
            return (D_pred - D_pre[i])^2
        end

        result = optimize(obj, -1000.0, 0.0)  # w bounds
        w_opt[i] = Optim.minimizer(result)
    end

    return w_opt
end
w_opt = pointwise_w_fit(vec(day_year), vec(EF), vec(POC3400), -3400.0, -100.0, 0.84)

P = compute_D(vec(day_year), vec(EF), -3400.0, -100.0, 0.84, w_opt)
#=
function compute_D0(t::Vector, D::Vector, z::Float64, z0::Float64, b::Float64, w::Float64)
    t_grid = reshape(t, 1, :)  # (1, Nt)
    Martin_term = ((z + z0) / (2*z0))^(-b)
    F_shifted = D ./ Martin_term
    # Compute the time before particles arrived at depth z
    t_unshifted = clamp.(t_grid .+ (z - z0) / w, minimum(t), maximum(t))
    # Interpolate back to D₀
    interp_F = linear_interpolation(t, F_shifted[:])  # flatten in case it's 2D
    D0_reconstructed = interp_F.(t_unshifted)
    return D0_reconstructed[:]
end
P0 = compute_D0(vec(day_year), vec(POC3400), -3400.0, -100.0, 0.84,-234.0)
=#
figure = Figure(size=(1000, 1000))

ax = Axis(figure[1, 1], title="Export flux 100 m", xlabel="Days after 1/1/2016", ylabel="flux (mg C m⁻² d⁻¹)")
lines!(ax, vec(day_year), vec(EF), color=:black, label = "True export flux")
# lines!(ax, vec(day_year), vec(P0), color=:red3, linestyle=:dash, label = "Reconstructed export flux")
axislegend(ax, position = (0.5, 1.0))

ax2 = Axis(figure[1, 2], title="POC flux 3400 m", xlabel="Days after 1/1/2016", ylabel="flux (mg C m⁻² d⁻¹)")
lines!(ax2, vec(day_year), vec(POC3400), color=:red3, label = "True POC flux 3400 m")
lines!(ax2, vec(day_year), vec(P), color=:black, linestyle=:dash, label = "Reconstructed flux at 3400 m")
axislegend(ax2, position = (0.7, 1.0))

ax3 = Axis(figure[1, 3], title="Optimized sinking speed (m d⁻¹)", xlabel="Days after 1/1/2016", ylabel="w (m d⁻¹)")
lines!(ax3, vec(day_year), vec(w_opt), color=:royalblue3, label = "sinking speed")

display(figure)

