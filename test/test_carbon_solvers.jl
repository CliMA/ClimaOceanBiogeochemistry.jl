using Test

include("../src/CarbonSystemSolvers.jl")

using .CarbonSystemSolvers
using .CarbonSystemSolvers:
    CarbonSystem,
    CarbonSystemParameters,
    CarbonSolverOptions,
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
        (log(coefficients.Cᵈⁱᶜₖ₁ᵣ₉₃),                               -13.4847, 4),
        (log(coefficients.Cᵈⁱᶜₖ₂ᵣ₉₃),                               -20.5504, 4),
        (log10(coefficients.Cᵈⁱᶜₖ₁ₗ₀₀),                              -5.8472, 4),
        (log10(coefficients.Cᵈⁱᶜₖ₂ₗ₀₀),                              -8.9660, 4),
        (log(coefficients.Cᵇₖ₁),                                    -19.7964, 4),
        (log(coefficients.Cᴴ²ᴼₖ₁ * coefficients.H⁺ₛoverH⁺ₜ) - 0.015, -30.434, 3),
        (log(coefficients.Cᴾᴼ⁴ₖ₁ * coefficients.H⁺ₛoverH⁺ₜ) - 0.015, -3.71, 2),
        (log(coefficients.Cᴾᴼ⁴ₖ₂ * coefficients.H⁺ₛoverH⁺ₜ) - 0.015, -13.727, 3),
        (log(coefficients.Cᴾᴼ⁴ₖ₃ * coefficients.H⁺ₛoverH⁺ₜ) - 0.015, -20.24, 2),
        (log(coefficients.Cˢⁱᵗₖ₁ * coefficients.H⁺ₛoverH⁺ₜ) - 0.015, -21.61, 2),
        (-log10(coefficients.Cᴴ²ˢₖ₁ * coefficients.H⁺ₛoverH⁺ₜ),       6.51, 2),
        (-log10(coefficients.Cᴺᴴ⁴ₖ₁ * coefficients.H⁺ₛoverH⁺ₜ),       9.26, 2),
        (log(coefficients.Cᴴᶠₖ₁),                                   -6.09, 2),
        (log(coefficients.Cᴴˢᴼ⁴ₖ₁),                                 -2.30, 2),
    ]

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

    # Regression values.
    #
    # These should ideally be generated once from an independently validated
    # configuration of CO2SYS, mocsy, or another documented reference solver.
    @test direct.pH ≈ 8.044006579188588 atol=pH_regression_atol rtol=0
    @test direct.pCO₂ᵒᶜᵉ * 1e6 ≈ 407.5213576262527 atol=pCO₂_regression_atol_μatm rtol=0

    @test correction.pH ≈ 8.036732866389743 atol=pH_regression_atol rtol=0
    @test correction.pCO₂ᵒᶜᵉ * 1e6 ≈ 415.09737958192645 atol=pCO₂_regression_atol_μatm rtol=0
              
    # NB, you wouldn't expect to get exactly the same results for UniversalRobustCarbonSystem 
    # because of the extra terms in the calcite alkalinity.
    @test universal.pH ≈ 8.037606899889317 atol=pH_regression_atol rtol=0
    @test universal.pCO₂ᵒᶜᵉ * 1e6 ≈ 414.18024641322165 atol=pCO₂_regression_atol_μatm rtol=0
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
# Representative environmental cases
# ------------------------------------------------------------------------

@testset "Universal solver across representative ocean conditions" begin
    carbon_params = CarbonSystemParameters()

    cases = [
        # Cold, high-latitude surface water.
        (;
            Θᶜ = -1.0,
            Sᴬ = 34.5,
            Δpᵦₐᵣ = 0.0,
            Cᵀ = 2150e-6,
            Aᵀ = 2300e-6,
            Pᵀ = 1.5e-6,
            Siᵀ = 10e-6,
        ),

        # Standard temperate surface ocean.
        (;
            Θᶜ = 25.0,
            Sᴬ = 35.0,
            Δpᵦₐᵣ = 0.0,
            Cᵀ = 2050e-6,
            Aᵀ = 2350e-6,
            Pᵀ = 0.5e-6,
            Siᵀ = 7.5e-6,
        ),

        # Warm, saline surface water.
        (;
            Θᶜ = 30.0,
            Sᴬ = 36.5,
            Δpᵦₐᵣ = 0.0,
            Cᵀ = 1950e-6,
            Aᵀ = 2400e-6,
            Pᵀ = 0.3e-6,
            Siᵀ = 2e-6,
        ),

        # Cold, nutrient-rich deep water at roughly 4000 dbar.
        (;
            Θᶜ = 2.0,
            Sᴬ = 35.0,
            Δpᵦₐᵣ = 400.0,
            Cᵀ = 2300e-6,
            Aᵀ = 2450e-6,
            Pᵀ = 2.0e-6,
            Siᵀ = 50e-6,
        ),
    ]

    for overrides in cases
        case = merge(nutrient_case, overrides)

        coefficients = CarbonChemistryCoefficients(
            carbon_params,
            case.Θᶜ,
            case.Sᴬ,
            case.Δpᵦₐᵣ,
        )

        solution = UniversalRobustCarbonSystem(; case...)

        assert_solution_invariants(solution, coefficients, case)

        # A deliberately broad plausibility bound for these ocean cases.
        @test 6.0 < solution.pH < 9.5
    end
end