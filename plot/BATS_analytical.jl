using GLMakie
using Printf
using Statistics
using Interpolations
using Optim
using XLSX

xlsx_data = XLSX.readxlsx("bats_flux.xlsx")
OC = xlsx_data["reconstruct"]

dayofyear = OC[2:end, 3]  # after 2020
OCflux_150 = OC[2:end, 4]
OCflux_300 = OC[2:end, 5]

function compute_D(t::Vector, D0::Vector, z::Float64, z0::Float64, b::Float64, w::Vector)
    t_grid = reshape(t, 1, :)   # (1, Nt)
    w_grid = reshape(w, 1, :)
    t_shifted = clamp.(t_grid .- (z .- z0) ./ w_grid, minimum(t), maximum(t))
    interp_D0 = linear_interpolation(t, D0)
    F_shifted = interp_D0.(t_shifted)
    Martin_term = ((z .+ z0) ./ (2*z0)).^(-b)
    return F_shifted .* Martin_term
end
#=
function objective(w_vec::Vector, t::Vector, D0::Vector, D_pre::Vector, z::Float64, z0::Float64, b::Float64)
    D_pred = compute_D(t, D0, z, z0, b, w_vec)
    return sum((D_pred .- D_pre).^2)
end

w_init = fill(-200.0, length(dayofyear))
# Define bounds 
lower = fill(-800.0, length(dayofyear))
upper = fill(0.0, length(dayofyear))

result = optimize(w -> objective(w, vec(dayofyear), vec(OCflux_150), vec(OCflux_300), -300.0, -150.0, 0.84),
                  lower, upper, w_init,
                  Fminbox(BFGS()))

w_opt = Optim.minimizer(result)
=#

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
w_opt = pointwise_w_fit(vec(dayofyear), vec(OCflux_150), vec(OCflux_300), -300.0, -150.0, 0.84)

P_300 = compute_D(vec(dayofyear), vec(OCflux_150), -300.0, -150.0, 0.84, w_opt)
# w_init[10:12] .= -1.0
# P_300_2 = compute_D(vec(dayofyear), vec(OCflux_150), -300.0, -150.0, 0.84, w_init)

############################################################
############# Reconstruct from deep to surface #############
############################################################

# function compute_D0(t::Vector, D::Vector, z::Float64, z0::Float64, b::Float64, w::Vector)
#     t_grid = reshape(t, 1, :)  # (1, Nt)
#     w_grid = reshape(w, 1, :)
#     Martin_term = ((z + z0) / (2*z0))^(-b)
#     F_shifted = D ./ Martin_term
#     # Compute the time before particles arrived at depth z
#     t_unshifted = clamp.(t_grid .+ (z - z0) ./ w_grid, minimum(t), maximum(t))
#     # Interpolate back to D₀
#     interp_F = linear_interpolation(t, F_shifted[:])  # flatten in case it's 2D
#     D0_reconstructed = interp_F.(t_unshifted)
#     return D0_reconstructed[:]
# end
# P_150 = compute_D0(vec(dayofyear), vec(OCflux_300), -300.0, -150.0, 0.84, wₛ)
#
figure = Figure(size=(1000, 1000))

# ax = Axis(figure[1, 1], title="Export flux 150 m", xlabel="Days after 1/1/2020", ylabel="flux (mg C m⁻² d⁻¹)")
# lines!(ax, vec(dayofyear), vec(OCflux_150), color=:black, label = "True export flux 150 m")
# # lines!(ax, vec(dayofyear), vec(P_150), color=:black, linestyle=:dash, label = "Reconstructed export flux 150 m")
# axislegend(ax, position = :rt)

ax2 = Axis(figure[1, 1], title="POC flux 300 m", xlabel="Days after 1/1/2020", ylabel="flux (mg C m⁻² d⁻¹)")
lines!(ax2, vec(dayofyear), vec(OCflux_300), color=:red3, label = "True POC flux 300 m")
lines!(ax2, vec(dayofyear), vec(P_300), color=:red3, linestyle=:dash, linewidth = 2,label = "Reconstructed flux")
# lines!(ax2, vec(dayofyear), vec(P_300_2), color=:blue, linestyle=:dot, linewidth = 2, label = "Reconstructed & adjusted w")
axislegend(ax2, position = :lt)
# xlims!(ax2, 300,500)

ax3 = Axis(figure[1, 2], title="Optimized sinking speed (m/d)", xlabel="Days after 1/1/2020", ylabel="w (m d⁻¹)")
lines!(ax3, vec(dayofyear), vec(w_opt), color=:royalblue3, label = "sinking speed")

display(figure)

