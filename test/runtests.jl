using Test
using ClimaOceanBiogeochemistry: CarbonSystemSolver

CarbonChemistryCoefficients() isa CarbonChemistryCoefficients

Θᶜ      = 25.0
Sᴬ      = 35.0
Δpᵦₐᵣ   = 0.0
Cᶜᵒᵉᶠᶠ = CarbonChemistryCoefficients(Θᶜ, Sᴬ, Δpᵦₐᵣ)

@assert round(log(Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀), digits = 4)     == -3.5617 # Handbook (2007)
@assert round(log(Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ᵣ₉₃), digits = 4)  == -13.4847 # Handbook (1994)
@assert round(log(Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ᵣ₉₃), digits = 4)  == -20.5504 # Handbook (1994)
#@assert Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ₘ₉₅  ==
#@assert Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ₘ₉₅  ==
@assert round(log10(Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ₗ₀₀), digits = 4)   == -5.8472 # Handbook (2007)
@assert round(log10(Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ₗ₀₀), digits = 4)   == -8.9660 # Handbook (2007)
@assert round(log(Cᶜᵒᵉᶠᶠ.Cᵇₖ₁), digits = 4)       == -19.7964 # Handbook (2007)
@assert round(log(Cᶜᵒᵉᶠᶠ.Cᴴ²ᴼₖ₁*Cᶜᵒᵉᶠᶠ.H⁺ₛoverH⁺ₜ)-0.015, digits = 3)     == -30.434 # Handbook (2007)
@assert round(log(Cᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁*Cᶜᵒᵉᶠᶠ.H⁺ₛoverH⁺ₜ)-0.015, digits = 2)     == -3.71 # Handbook (2007)
@assert round(log(Cᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂*Cᶜᵒᵉᶠᶠ.H⁺ₛoverH⁺ₜ)-0.015, digits = 3)     == -13.727 # Handbook (2007)
@assert round(log(Cᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₃*Cᶜᵒᵉᶠᶠ.H⁺ₛoverH⁺ₜ)-0.015, digits = 2)     == -20.24 # Handbook (2007)
@assert round(log(Cᶜᵒᵉᶠᶠ.Cˢⁱᵗₖ₁*Cᶜᵒᵉᶠᶠ.H⁺ₛoverH⁺ₜ)-0.015, digits = 2)     == -21.61 # Handbook (2007)
@assert round(-log10(Cᶜᵒᵉᶠᶠ.Cᴴ²ˢₖ₁*Cᶜᵒᵉᶠᶠ.H⁺ₛoverH⁺ₜ), digits = 2)     == 6.51 # Lewis and Wallace (1998)
@assert round(-log10(Cᶜᵒᵉᶠᶠ.Cᴺᴴ⁴ₖ₁*Cᶜᵒᵉᶠᶠ.H⁺ₛoverH⁺ₜ), digits = 2)     == 9.26 # Lewis and Wallace (1998)
#@assert Cᶜᵒᵉᶠᶠ.Cᴴᶠᵦ₁      ==
@assert round(log(Cᶜᵒᵉᶠᶠ.Cᴴᶠₖ₁), digits = 2)       == -6.09 # Handbook (2007)
@assert round(log(Cᶜᵒᵉᶠᶠ.Cᴴˢᴼ⁴ₖ₁), digits = 2)     == -2.30 # Handbook (2007)
#@assert Cᶜᵒᵉᶠᶠ.Cᶜᵃˡᶜⁱᵗᵉₛₚ  ==
#@assert Cᶜᵒᵉᶠᶠ.Cᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ==
#@assert Cᶜᵒᵉᶠᶠ.Cᴮᵀ        ==
#@assert Cᶜᵒᵉᶠᶠ.Cᶠᵀ        ==
#@assert Cᶜᵒᵉᶠᶠ.Cᶜᵃ        ==
#@assert Cᶜᵒᵉᶠᶠ.Cˢᴼ⁴       ==