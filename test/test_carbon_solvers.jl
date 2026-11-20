using Test
#using Printf

include("../src/CarbonSystemSolvers.jl")

using .CarbonSystemSolvers
using .CarbonSystemSolvers:
    CarbonSystem,
    CarbonSystemParameters,
    CarbonSolverOptions,
    CarbonCalculationConstants,
    CarbonCoefficientParameters,
    CarbonChemistryCoefficients

using .CarbonSystemSolvers.DirectCubicCarbonSolver:
    DirectCubicCarbonSystem

using .CarbonSystemSolvers.AlkalinityCorrectionCarbonSolver:
    AlkalinityCorrectionCarbonSystem

using .CarbonSystemSolvers.UniversalRobustCarbonSolver:
    UniversalRobustCarbonSystem

# ------------------------------------------------------------------------
# Test parameters
# ------------------------------------------------------------------------

Θᶜ = 25.0
Sᴬ = 35.0
Δpᵦₐᵣ = 0.0       # applied pressure [bar]
Cᵀ = 2050e-6      # [mol kg⁻¹]
Aᵀ = 2350e-6      # [mol kg⁻¹]
Pᵀ = 0.5e-6       # [mol kg⁻¹]
Siᵀ = 7.5e-6      # [mol kg⁻¹]
pH₀ = 8.0          # Initial solver guess; do not overwrite this.
pCO₂ᵃᵗᵐ = 280e-6  # [atm]

base_case = (;
    Θᶜ,
    Sᴬ,
    Δpᵦₐᵣ,
    Cᵀ,
    Aᵀ,
    pH = pH₀,
    pCO₂ᵃᵗᵐ,
)

nutrient_case = merge(base_case, (; Pᵀ, Siᵀ))

# These are deliberately absolute tolerances in physically meaningful units.
#
# They are much tighter than scientifically meaningful uncertainty, but allow
# harmless differences from compiler, architecture, or operation ordering.
pH_regression_atol = 1e-8
pCO₂_regression_atol_μatm = 1e-5

# Closure tolerances should be near floating-point roundoff, rather than the
# much looser regression tolerances above.
carbon_closure_rtol = 5e-12
carbon_closure_atol = 1e-15
henry_closure_rtol = 5e-12
henry_closure_atol = 1e-16

# ------------------------------------------------------------------------
# Parameter and coefficient tests
# ------------------------------------------------------------------------

@testset "Carbon-system parameters" begin
    carbon_params = CarbonSystemParameters()

    @test carbon_params isa CarbonSystemParameters
    @test carbon_params.Sᵒᵖᵗˢ isa CarbonSolverOptions
    @test carbon_params.Pᵈⁱᶜₖ₀ isa CarbonCoefficientParameters
end

@testset "Carbon chemistry coefficients" begin
    carbon_params = CarbonSystemParameters()
    coefficients = CarbonChemistryCoefficients(carbon_params, Θᶜ, Sᴬ, Δpᵦₐᵣ)

    @test coefficients isa CarbonChemistryCoefficients

    # The references below are reported to a finite number of decimal places.
    # Test against half a final decimal place, rather than using `round(x) == y`.
    #
    # Each tuple is: (computed value, literature value, decimal places reported).
    coefficient_checks = [
        (log(coefficients.Cᵈⁱᶜₖ₀),                                  -3.5617, 4),
        (log(coefficients.Cᵇₖ₁),                                    -19.7964, 4),
        (log(coefficients.Cᴴ²ᴼₖ₁ * coefficients.H⁺ₛoverH⁺ₜ) - 0.015, -30.434, 3),
        (log(coefficients.Cᴾᴼ⁴ₖ₁ * coefficients.H⁺ₛoverH⁺ₜ) - 0.015, -3.71, 2),
        (log(coefficients.Cᴾᴼ⁴ₖ₂ * coefficients.H⁺ₛoverH⁺ₜ) - 0.015, -13.727, 3),
        (log(coefficients.Cᴾᴼ⁴ₖ₃ * coefficients.H⁺ₛoverH⁺ₜ) - 0.015, -20.24, 2),
        (log(coefficients.Cˢⁱᵗₖ₁ * coefficients.H⁺ₛoverH⁺ₜ) - 0.015, -21.61, 2),
        (-log10(coefficients.Cᴴ²ˢₖ₁ * coefficients.H⁺ₛoverH⁺ₜ),       6.51, 2),
        (-log10(coefficients.Cᴺᴴ⁴ₖ₁ * coefficients.H⁺ₛoverH⁺ₜ),       9.26, 2),
        (log(coefficients.Cᴴˢᴼ⁴ₖ₁),                                 -2.30, 2),
    ]

    # Calculate both Roy1993 and Lueker2000 carbonic acid coefficients for comparison
    carbon_params = CarbonSystemParameters(; carbonic_acid = Roy1993())
    coefficients = CarbonChemistryCoefficients(carbon_params, Θᶜ, Sᴬ, Δpᵦₐᵣ)

    push!(coefficient_checks,
        (log(coefficients.Cᵈⁱᶜₖ₁),                               -13.4847, 4),
        (log(coefficients.Cᵈⁱᶜₖ₂),                               -20.5504, 4),
    )

    carbon_params = CarbonSystemParameters(; carbonic_acid = Lueker2000())
    coefficients = CarbonChemistryCoefficients(carbon_params, Θᶜ, Sᴬ, Δpᵦₐᵣ)

    push!(coefficient_checks,
        (log10(coefficients.Cᵈⁱᶜₖ₁),                              -5.8472, 4),
        (log10(coefficients.Cᵈⁱᶜₖ₂),                              -8.9660, 4),
    )

    carbon_params = CarbonSystemParameters(; fluoric_acid = PerezFraga1987())
    coefficients = CarbonChemistryCoefficients(carbon_params, Θᶜ, Sᴬ, Δpᵦₐᵣ)

    push!(coefficient_checks,
            (log(coefficients.Cᴴᶠₖ₁),                             -6.09, 2),
        )

    for (actual, reference, digits) in coefficient_checks
        reference_atol = 0.5001 * 10.0^(-digits)
        @test actual ≈ reference atol=reference_atol rtol=0
    end
end

# ------------------------------------------------------------------------
# Constructor and standard-reference tests
# ------------------------------------------------------------------------

@testset "Carbon-system constructors" begin
    @test DirectCubicCarbonSystem() isa CarbonSystem
    @test AlkalinityCorrectionCarbonSystem() isa CarbonSystem
    @test UniversalRobustCarbonSystem() isa CarbonSystem
end

"""
Check invariants that every returned carbonate-system solution should satisfy.

The K₀ relation assumes that `pCO₂ᵒᶜᵉ` is the partial-pressure variable
associated with the solver's `Cᵈⁱᶜₖ₀` convention.
"""
function assert_solution_invariants(solution, coefficients, inputs)

    finite_values = (
        solution.pH,
        solution.CO₂ˢᵒˡ,
        solution.HCO₃⁻,
        solution.CO₃²⁻,
        solution.Cᵀ,
        solution.Aᵀ,
        solution.pCO₂ᵒᶜᵉ,
        solution.pCO₂ᵃᵗᵐ,
    )

    @test all(isfinite, finite_values)

    positive_values = (
        solution.CO₂ˢᵒˡ,
        solution.HCO₃⁻,
        solution.CO₃²⁻,
        solution.Cᵀ,
        solution.Aᵀ,
        solution.pCO₂ᵒᶜᵉ,
        solution.pCO₂ᵃᵗᵐ,
    )

    @test all(x -> x > 0, positive_values)
    @test solution.pH > 0

    # Inputs should be preserved in the returned solution.
    @test solution.Cᵀ ≈ inputs.Cᵀ rtol=carbon_closure_rtol atol=carbon_closure_atol
    @test solution.Aᵀ ≈ inputs.Aᵀ rtol=carbon_closure_rtol atol=carbon_closure_atol
    @test solution.pCO₂ᵃᵗᵐ == inputs.pCO₂ᵃᵗᵐ

    # DIC partitioning must close.
    carbon_sum = solution.CO₂ˢᵒˡ + solution.HCO₃⁻ + solution.CO₃²⁻

    @test carbon_sum ≈ solution.Cᵀ rtol=carbon_closure_rtol atol=carbon_closure_atol

    # Henry-law closure.
    @test solution.CO₂ˢᵒˡ ≈ coefficients.Cᵈⁱᶜₖ₀ * solution.pCO₂ᵒᶜᵉ rtol=henry_closure_rtol atol=henry_closure_atol

    return nothing
end

# AlkalinityCorrectionCarbonSystem should be iterated to convergence from a cold 
# start, but is accurate enough for a warm start from previous timestep pH values.
function iterated_alkalinity_correction(;
        Θᶜ,
        Sᴬ,
        Δpᵦₐᵣ,
        Cᵀ,
        Aᵀ,
        Pᵀ,
        Siᵀ,
        pH,
        pCO₂ᵃᵗᵐ,
        pH_atol = 1e-10,
        maxiter = 100,
    )

    pH_guess = pH

    for iteration in 1:maxiter
        solution = AlkalinityCorrectionCarbonSystem(
                Θᶜ      = Θᶜ,
                Sᴬ      = Sᴬ,
                Δpᵦₐᵣ   = Δpᵦₐᵣ,
                Cᵀ      = Cᵀ,
                Aᵀ      = Aᵀ,
                Pᵀ      = Pᵀ,
                Siᵀ     = Siᵀ,
                pH      = pH_guess,
                pCO₂ᵃᵗᵐ = pCO₂ᵃᵗᵐ,
        )

        pH_new = solution.pH
        residual = pH_new - pH_guess

        if abs(residual) ≤ pH_atol
            return solution, iteration, residual
        end

        pH_guess = pH_new
    end

    error(
        "AlkalinityCorrectionCarbonSystem failed to converge after $maxiter iterations. " *
        "Initial pH = $pH₀, final pH iterate = $pH_guess."
    )
end

@testset "Standard carbonate-system reference case" begin
    carbon_params = CarbonSystemParameters()
    coefficients = CarbonChemistryCoefficients(carbon_params, Θᶜ, Sᴬ, Δpᵦₐᵣ)

    direct = DirectCubicCarbonSystem(; base_case...)
    universal = UniversalRobustCarbonSystem(; nutrient_case...)
    correction, iteration, residual = iterated_alkalinity_correction(; nutrient_case...)

    # Structural and thermodynamic consistency tests.
    assert_solution_invariants(direct, coefficients, base_case)
    assert_solution_invariants(correction, coefficients, nutrient_case)
    assert_solution_invariants(universal, coefficients, nutrient_case)

    # Regression values. You wouldn't expect to get exactly the same results 
    # between solvers because of the extra terms in the calcite alkalinity.
    @test direct.pH ≈ 8.044006579188588 atol=pH_regression_atol rtol=0
    @test direct.pCO₂ᵒᶜᵉ * 1e6 ≈ 407.5213576262527 atol=pCO₂_regression_atol_μatm rtol=0

    @test correction.pH ≈ 8.036732847715635 atol=pH_regression_atol rtol=0                   
    @test correction.pCO₂ᵒᶜᵉ * 1e6 ≈ 415.097399196748 atol=pCO₂_regression_atol_μatm rtol=0
    
    @test universal.pH ≈ 8.036732758889796 atol=pH_regression_atol rtol=0
    @test universal.pCO₂ᵒᶜᵉ * 1e6 ≈ 415.09749249722745 atol=pCO₂_regression_atol_μatm rtol=0
end

# ------------------------------------------------------------------------
# Solver-behaviour tests
# ------------------------------------------------------------------------

@testset "Universal solver is insensitive to initial pH guess" begin
    reference = UniversalRobustCarbonSystem(; nutrient_case...)

    for initial_pH in (7.0, 8.0, 9.0)
        case = merge(nutrient_case, (; pH = initial_pH))
        solution = UniversalRobustCarbonSystem(; case...)

        @test solution.pH ≈ reference.pH atol=pH_regression_atol rtol=0

        @test solution.pCO₂ᵒᶜᵉ * 1e6 ≈ reference.pCO₂ᵒᶜᵉ * 1e6 atol=pCO₂_regression_atol_μatm rtol=0
    end
end

@testset "Atmospheric pCO₂ does not alter closed-system speciation" begin
    low_air_case = merge(nutrient_case, (; pCO₂ᵃᵗᵐ = 180e-6))
    high_air_case = merge(nutrient_case, (; pCO₂ᵃᵗᵐ = 600e-6))

    low_air = UniversalRobustCarbonSystem(; low_air_case...)
    high_air = UniversalRobustCarbonSystem(; high_air_case...)

    @test low_air.pCO₂ᵃᵗᵐ == 180e-6
    @test high_air.pCO₂ᵃᵗᵐ == 600e-6

    @test low_air.pH ≈ high_air.pH atol=pH_regression_atol rtol=0

    @test low_air.pCO₂ᵒᶜᵉ * 1e6 ≈ high_air.pCO₂ᵒᶜᵉ * 1e6 atol=pCO₂_regression_atol_μatm rtol=0
end

@testset "Expected carbonate-system sensitivities" begin
    # At fixed alkalinity, adding DIC should lower pH and raise pCO₂.
    low_DIC_case = merge(nutrient_case, (; Cᵀ = 1950e-6))
    high_DIC_case = merge(nutrient_case, (; Cᵀ = 2150e-6))

    low_DIC = UniversalRobustCarbonSystem(; low_DIC_case...)
    high_DIC = UniversalRobustCarbonSystem(; high_DIC_case...)

    @test high_DIC.pH < low_DIC.pH
    @test high_DIC.pCO₂ᵒᶜᵉ > low_DIC.pCO₂ᵒᶜᵉ

    # At fixed DIC, adding alkalinity should raise pH and lower pCO₂.
    low_alk_case = merge(nutrient_case, (; Aᵀ = 2250e-6))
    high_alk_case = merge(nutrient_case, (; Aᵀ = 2450e-6))

    low_alk = UniversalRobustCarbonSystem(; low_alk_case...)
    high_alk = UniversalRobustCarbonSystem(; high_alk_case...)

    @test high_alk.pH > low_alk.pH
    @test high_alk.pCO₂ᵒᶜᵉ < low_alk.pCO₂ᵒᶜᵉ
end

# ------------------------------------------------------------------------
# Test derivative consistency of the UniversalRobustCarbonSolver
# ------------------------------------------------------------------------
@testset "Derivative consistency of the UniversalRobustCarbonSolver" begin
    carbon_params = CarbonSystemParameters(
        carbonic_acid = Millero1995()
    )

    coefficients = CarbonChemistryCoefficients(
        carbon_params, 25.0, 35.0, 0.0
    )

    coeff = (
        Cᴾᴼ⁴ₖ₁ = coefficients.Cᴾᴼ⁴ₖ₁,
        Cᴾᴼ⁴ₖ₂ = coefficients.Cᴾᴼ⁴ₖ₂,
        Cᴾᴼ⁴ₖ₃ = coefficients.Cᴾᴼ⁴ₖ₃,
    )

    Pᵀ = 0.5e-6

    for pH in (5.3, 7.0, 8.0, 9.0, 11.0)
        H = 10.0^-pH
        δH = H * 1e-5

        numerical = (
            CarbonSystemSolvers.UniversalRobustCarbonSolver.FAPᵀ(Pᵀ, H + δH, coeff) -
            CarbonSystemSolvers.UniversalRobustCarbonSolver.FAPᵀ(Pᵀ, H - δH, coeff)
        ) / (2δH)

        analytical = CarbonSystemSolvers.UniversalRobustCarbonSolver.F∂A∂Pᵀ(Pᵀ, H, coeff)

        @test isapprox(
            numerical, analytical;
            rtol = 1e-7,
            atol = 1e-12
        )
    end
end

# ------------------------------------------------------------------------
# Test UniversalRobustCarbonSolver with different (and extreme) environmental cases
# ------------------------------------------------------------------------
@testset "Original Fortran SolveSAPHE reference cases" begin
    cases = [
        (; id = 1, label = "Warm surface - baseline", source = "atgen_sw1-sw.bin",
           Θᶜ = 25.0, Sᴬ = 35.0, Δpᵦₐᵣ = 0.0, pCO₂ᵃᵗᵐ = 280.0e-6,
           Cᵀ = 2050.5e-6, Aᵀ = 2350.5e-6, Pᵀ = 0.5e-6, Siᵀ = 5.0e-6,
           pH_SWS = 7.9996792918870065, H_SWS = 1.00073873044721231e-08),
        (; id = 2, label = "Warm surface - high DIC", source = "atgen_sw1-sw.bin",
           Θᶜ = 25.0, Sᴬ = 35.0, Δpᵦₐᵣ = 0.0, pCO₂ᵃᵗᵐ = 280.0e-6,
           Cᵀ = 2400.5e-6, Aᵀ = 2300.5e-6, Pᵀ = 0.5e-6, Siᵀ = 5.0e-6,
           pH_SWS = 7.0414258788317188, H_SWS = 9.09021428784569179e-08),
        (; id = 3, label = "Cold surface - baseline", source = "atgen_sw1-sc.bin",
           Θᶜ = 2.0, Sᴬ = 35.0, Δpᵦₐᵣ = 0.0, pCO₂ᵃᵗᵐ = 280.0e-6,
           Cᵀ = 2150.5e-6, Aᵀ = 2300.5e-6, Pᵀ = 0.5e-6, Siᵀ = 5.0e-6,
           pH_SWS = 8.0808809806938697, H_SWS = 8.30078221009040016e-09),
        (; id = 4, label = "Cold surface - high alkalinity", source = "atgen_sw1-sc.bin",
           Θᶜ = 2.0, Sᴬ = 35.0, Δpᵦₐᵣ = 0.0, pCO₂ᵃᵗᵐ = 280.0e-6,
           Cᵀ = 2050.5e-6, Aᵀ = 2450.5e-6, Pᵀ = 0.5e-6, Siᵀ = 5.0e-6,
           pH_SWS = 8.5294619536674947, H_SWS = 2.95486773808702780e-09),
        (; id = 5, label = "Deep cold - matched baseline", source = "atgen_sw1-dc.bin",
           Θᶜ = 2.0, Sᴬ = 35.0, Δpᵦₐᵣ = 300.0, pCO₂ᵃᵗᵐ = 280.0e-6,
           Cᵀ = 2150.5e-6, Aᵀ = 2300.5e-6, Pᵀ = 0.5e-6, Siᵀ = 5.0e-6,
           pH_SWS = 7.9633874789520078, H_SWS = 1.08795897963544235e-08),
        (; id = 6, label = "Deep cold - matched high alkalinity", source = "atgen_sw1-dc.bin",
           Θᶜ = 2.0, Sᴬ = 35.0, Δpᵦₐᵣ = 300.0, pCO₂ᵃᵗᵐ = 280.0e-6,
           Cᵀ = 2050.5e-6, Aᵀ = 2450.5e-6, Pᵀ = 0.5e-6, Siᵀ = 5.0e-6,
           pH_SWS = 8.4181475037424693, H_SWS = 3.81814569346148547e-09),
        (; id = 7, label = "Brackish cold - matched baseline", source = "atgen_sw1-sb.bin",
           Θᶜ = 2.0, Sᴬ = 3.5, Δpᵦₐᵣ = 0.0, pCO₂ᵃᵗᵐ = 280.0e-6,
           Cᵀ = 2150.5e-6, Aᵀ = 2300.5e-6, Pᵀ = 0.5e-6, Siᵀ = 5.0e-6,
           pH_SWS = 8.6649964661595593, H_SWS = 2.16273612177073106e-09),
        (; id = 8, label = "Brackish cold - matched high alkalinity", source = "atgen_sw1-sb.bin",
           Θᶜ = 2.0, Sᴬ = 3.5, Δpᵦₐᵣ = 0.0, pCO₂ᵃᵗᵐ = 280.0e-6,
           Cᵀ = 2050.5e-6, Aᵀ = 2450.5e-6, Pᵀ = 0.5e-6, Siᵀ = 5.0e-6,
           pH_SWS = 9.1525857159523589, H_SWS = 7.03743317624632166e-10),
        (; id = 9, label = "Warm extreme - high DIC/low alkalinity", source = "atgen_sw3-sw.bin",
           Θᶜ = 25.0, Sᴬ = 35.0, Δpᵦₐᵣ = 0.0, pCO₂ᵃᵗᵐ = 280.0e-6,
           Cᵀ = 4505.0e-6, Aᵀ = 1005.0e-6, Pᵀ = 0.5e-6, Siᵀ = 5.0e-6,
           pH_SWS = 5.3071907421209330, H_SWS = 4.92957249596676450e-06),
        (; id = 10, label = "Warm extreme - low DIC/high alkalinity", source = "atgen_sw3-sw.bin",
           Θᶜ = 25.0, Sᴬ = 35.0, Δpᵦₐᵣ = 0.0, pCO₂ᵃᵗᵐ = 280.0e-6,
           Cᵀ = 1005.0e-6, Aᵀ = 4005.0e-6, Pᵀ = 0.5e-6, Siᵀ = 5.0e-6,
           pH_SWS = 10.4176074925104096, H_SWS = 3.82289621193899265e-11),
        (; id = 11, label = "Dilute cold - equal DIC/alkalinity", source = "atgen_sw4-sc.bin",
           Θᶜ = 2.0, Sᴬ = 35.0, Δpᵦₐᵣ = 0.0, pCO₂ᵃᵗᵐ = 280.0e-6,
           Cᵀ = 505.0e-6, Aᵀ = 505.0e-6, Pᵀ = 0.5e-6, Siᵀ = 5.0e-6,
           pH_SWS = 7.4349088936928442, H_SWS = 3.67359357070350254e-08),
        (; id = 12, label = "Dilute cold - low DIC/high alkalinity", source = "atgen_sw4-sc.bin",
           Θᶜ = 2.0, Sᴬ = 35.0, Δpᵦₐᵣ = 0.0, pCO₂ᵃᵗᵐ = 280.0e-6,
           Cᵀ = 105.0e-6, Aᵀ = 1405.0e-6, Pᵀ = 0.5e-6, Siᵀ = 5.0e-6,
           pH_SWS = 11.0938957608298576, H_SWS = 8.05571770978109956e-12),
    ]
    
    params = CarbonSystemParameters(; 
                                Cᶜᵒⁿˢᵗ = CarbonCalculationConstants(
                                        R = 83.14510,
                                ),
                                carbonic_acid = Millero1995(),
                                silicic_acid  = YaoMillero1995(),
    )

    #println("\nCase  Julia pH(T)   Julia pH(SWS) Fortran pH(SWS)      ΔpH(SWS)")
    for case in cases
        (; Θᶜ, Sᴬ, Δpᵦₐᵣ, Cᵀ, Aᵀ, Pᵀ, Siᵀ, pCO₂ᵃᵗᵐ) = case
        coefficients = CarbonChemistryCoefficients(params, Θᶜ, Sᴬ, Δpᵦₐᵣ)
        solution     = UniversalRobustCarbonSystem(; Θᶜ, Sᴬ, Δpᵦₐᵣ, Cᵀ, Aᵀ, Pᵀ, Siᵀ,
                                               pH = 8.0, pCO₂ᵃᵗᵐ = pCO₂ᵃᵗᵐ, params)

        assert_solution_invariants(solution, coefficients, case)

        # Reference SolveSAPHE output is on the seawater scale (SWS)
        # [H+]_SWS / [H+]_Total = coefficients.H⁺ₛoverH⁺ₜ
        # Therefore pH_SWS = pH_Total - log10(H_SWS/H_Total).
        julia_pH_SWS = solution.pH - log10(coefficients.H⁺ₛoverH⁺ₜ)

    #   δpH = julia_pH_SWS - case.pH_SWS
    #    @printf("%4d  %12.9f  %14.9f  %15.9f  %+13.7e  %s\n",
    #            case.id, solution.pH, julia_pH_SWS, case.pH_SWS, δpH, case.label)
        @test isfinite(julia_pH_SWS)

        @test julia_pH_SWS ≈ case.pH_SWS atol = 1e-12 rtol = 0
    end
end