using Test

include("../src/CarbonSystemSolvers.jl")

using .CarbonSystemSolvers
using .CarbonSystemSolvers: CarbonCoefficientParameters, 
                             CarbonSolverParameters, 
                             CarbonSystemParameters,
                             CarbonSystem, 
                             CarbonChemistryCoefficients
using .CarbonSystemSolvers.DirectCubicCarbonSolver: DirectCubicCarbonSystem
using .CarbonSystemSolvers.AlkalinityCorrectionCarbonSolver: AlkalinityCorrectionCarbonSystem
using .CarbonSystemSolvers.UniversalRobustCarbonSolver: UniversalRobustCarbonSystem

# Test carbon chemistry coefficients
Θᶜ      = 25.0
Sᴬ      = 35.0
Δpᵦₐᵣ   = 0.0
Cᵀ      = 2050e-6 # umol/kg to mol/kg
Aᵀ      = 2350e-6 # umol/kg to mol/kg
pCO₂ᵃᵗᵐ = 280e-6  # uatm to atm
pH      = 8.0

# Use default carbon system parameters
carbon_params = CarbonSystemParameters()

@test carbon_params       isa CarbonSystemParameters
#@test carbon_params.Sᵒᵖᵗˢ isa CarbonSolverParameters
#@test carbon_params.Pᵈⁱᶜₖ₀ isa CarbonCoefficientParameters

Cᶜᵒᵉᶠᶠ = CarbonChemistryCoefficients(carbon_params, Θᶜ, Sᴬ, Δpᵦₐᵣ)

# Check the type of the coefficients struct
@test Cᶜᵒᵉᶠᶠ isa CarbonChemistryCoefficients

#Check the values of calculated constants (sources in comments)
@test round(log(Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀),                        digits = 4) == -3.5617 # Handbook (2007)
@test round(log(Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ᵣ₉₃),                     digits = 4) == -13.4847 # Handbook (1994)
@test round(log(Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ᵣ₉₃),                     digits = 4) == -20.5504 # Handbook (1994)
#@test Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ₘ₉₅  ≈
#@test Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ₘ₉₅  ≈
@test round(log10(Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ₗ₀₀),                    digits = 4) == -5.8472 # Handbook (2007)
@test round(log10(Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ₗ₀₀),                    digits = 4) == -8.9660 # Handbook (2007)
@test round(log(Cᶜᵒᵉᶠᶠ.Cᵇₖ₁),                          digits = 4) == -19.7964 # Handbook (2007)
@test round(log(Cᶜᵒᵉᶠᶠ.Cᴴ²ᴼₖ₁*Cᶜᵒᵉᶠᶠ.H⁺ₛoverH⁺ₜ)-0.015, digits = 3) == -30.434 # Handbook (2007)
@test round(log(Cᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁*Cᶜᵒᵉᶠᶠ.H⁺ₛoverH⁺ₜ)-0.015, digits = 2) == -3.71 # Handbook (2007)
@test round(log(Cᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂*Cᶜᵒᵉᶠᶠ.H⁺ₛoverH⁺ₜ)-0.015, digits = 3) == -13.727 # Handbook (2007)
@test round(log(Cᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₃*Cᶜᵒᵉᶠᶠ.H⁺ₛoverH⁺ₜ)-0.015, digits = 2) == -20.24 # Handbook (2007)
@test round(log(Cᶜᵒᵉᶠᶠ.Cˢⁱᵗₖ₁*Cᶜᵒᵉᶠᶠ.H⁺ₛoverH⁺ₜ)-0.015, digits = 2) == -21.61 # Handbook (2007)
@test round(-log10(Cᶜᵒᵉᶠᶠ.Cᴴ²ˢₖ₁*Cᶜᵒᵉᶠᶠ.H⁺ₛoverH⁺ₜ),    digits = 2) == 6.51 # Lewis and Wallace (1998)
@test round(-log10(Cᶜᵒᵉᶠᶠ.Cᴺᴴ⁴ₖ₁*Cᶜᵒᵉᶠᶠ.H⁺ₛoverH⁺ₜ),    digits = 2) == 9.26 # Lewis and Wallace (1998)
#@test Cᶜᵒᵉᶠᶠ.Cᴴᶠᵦ₁      ≈
@test round(log(Cᶜᵒᵉᶠᶠ.Cᴴᶠₖ₁),                         digits = 2) == -6.09 # Handbook (2007)
@test round(log(Cᶜᵒᵉᶠᶠ.Cᴴˢᴼ⁴ₖ₁),                       digits = 2) == -2.30 # Handbook (2007)
#@test Cᶜᵒᵉᶠᶠ.Cᶜᵃˡᶜⁱᵗᵉₛₚ  ≈
#@test Cᶜᵒᵉᶠᶠ.Cᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ≈
#@test Cᶜᵒᵉᶠᶠ.Cᴮᵀ        ≈
#@test Cᶜᵒᵉᶠᶠ.Cᶠᵀ        ≈
#@test Cᶜᵒᵉᶠᶠ.Cᶜᵃ        ≈
#@test Cᶜᵒᵉᶠᶠ.Cˢᴼ⁴       ≈

# Test the CarbonSystemSolvers module
@test DirectCubicCarbonSystem() isa CarbonSystem

(; pH, CO₂ˢᵒˡ, HCO₃⁻, CO₃²⁻, Cᵀ, Aᵀ, pCO₂ᵒᶜᵉ, pCO₂ᵃᵗᵐ) = 
DirectCubicCarbonSystem(
        Θᶜ      = Θᶜ,
        Sᴬ      = Sᴬ,
        Δpᵦₐᵣ   = Δpᵦₐᵣ,
        Cᵀ      = Cᵀ,
        Aᵀ      = Aᵀ,
        pH      = pH,
        pCO₂ᵃᵗᵐ = pCO₂ᵃᵗᵐ,
        )

@test pH            == 8.044006579710093 
@test pCO₂ᵃᵗᵐ * 1e6 == 280.0             
@test pCO₂ᵒᶜᵉ * 1e6 == 407.52135708764496

Pᵀ = 0.5e-6  # umol/kg to mol/kg
Siᵀ = 7.5e-6 # umol/kg to mol/kg

@test AlkalinityCorrectionCarbonSystem() isa CarbonSystem

(; pH, CO₂ˢᵒˡ, HCO₃⁻, CO₃²⁻, Cᵀ, Aᵀ, pCO₂ᵒᶜᵉ, pCO₂ᵃᵗᵐ) = 
AlkalinityCorrectionCarbonSystem(
        Θᶜ      = Θᶜ,
        Sᴬ      = Sᴬ,
        Δpᵦₐᵣ   = Δpᵦₐᵣ,
        Cᵀ      = Cᵀ,
        Aᵀ      = Aᵀ,
        Pᵀ      = Pᵀ,
        Siᵀ     = Siᵀ,
        pH      = pH,
        pCO₂ᵃᵗᵐ = pCO₂ᵃᵗᵐ,
        )

@test pH            == 8.033988293659919
@test pCO₂ᵃᵗᵐ * 1e6 == 280.0            
@test pCO₂ᵒᶜᵉ * 1e6 == 417.9894057400246

@test UniversalRobustCarbonSystem() isa CarbonSystem

(; pH, CO₂ˢᵒˡ, HCO₃⁻, CO₃²⁻, Cᵀ, Aᵀ, pCO₂ᵒᶜᵉ, pCO₂ᵃᵗᵐ) = 
UniversalRobustCarbonSystem(
        Θᶜ      = Θᶜ,
        Sᴬ      = Sᴬ,
        Δpᵦₐᵣ   = Δpᵦₐᵣ,
        Cᵀ      = Cᵀ,
        Aᵀ      = Aᵀ,
        Pᵀ      = Pᵀ,
        Siᵀ     = Siᵀ,
        pH      = pH,
        pCO₂ᵃᵗᵐ = pCO₂ᵃᵗᵐ,
        )

@test pH            == 8.037606899889317
@test pCO₂ᵃᵗᵐ * 1e6 == 280.0            
@test pCO₂ᵒᶜᵉ * 1e6 == 414.18024641322165
# NB, you wouldn't expect to get exactly the same results here 
# because of the extra terms in the calcite alkalinity.
