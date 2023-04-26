module CarbonSystemSolvers
export CarbonSystem

struct CarbonSystem{FT}
    pH     :: FT
    CO₂ˢᵒˡ :: FT
    HCO₃⁻  :: FT
    CO₃²⁻  :: FT
    Cᵀ     :: FT
    Aᵀ     :: FT
    pCO₂ᵒᶜᵉ:: FT
    pCO₂ᵃᵗᵐ:: FT
end

"""
CarbonSolverApprox solves a cubic equation in terms of [H⁺]; 
Not for serious use, but as a placeholder and for testing purposes
"""
module CarbonSolverApprox
export CarbonSystemApprox

using RootSolvers
using ..CarbonSystemSolvers: CarbonSystem
include("carbon_chemistry_coefficients.jl")

"""
CarbonSystemApprox(
        Θ       :: FT = 25.0,
        Sᴬ      :: FT = 35.0,
        Δpᵦₐᵣ   :: FT = 0.0,
        Cᵀ      :: FT = 2050.0e-6,
        Aᵀ      :: FT = 2350.0e-6,
        pH      :: FT = 8.0,
        pCO₂ᵃᵗᵐ :: FT = 280.0e-6,
        )

TBW
"""
@inline function CarbonSystemApprox(
        Θᶜ      :: FT = 25.0,
        Sᴬ      :: FT = 35.0,
        Δpᵦₐᵣ   :: FT = 0.0,
        Cᵀ      :: FT = 2050.0e-6,
        Aᵀ      :: FT = 2350.0e-6,
        pH      :: FT = 8.0,
        pCO₂ᵃᵗᵐ :: FT = 280.0e-6) where {FT}

    # CarbonChemistryCoefficients are pretty much all in mol/kg, hence the 1e-6 factors for Cᵀ and Aᵀ
    Cᶜᵒᵉᶠᶠ = CarbonChemistryCoefficients(Θᶜ, Sᴬ, Δpᵦₐᵣ)
    
    # Some logic here about choosing coefficient options, particularly Cᵈⁱᶜ 
    Pᶜᵒᵉᶠᶠ = (Cᵈⁱᶜₖ₀ = Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀,
              Cᵈⁱᶜₖ₁ = Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ₗ₀₀,
              Cᵈⁱᶜₖ₂ = Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ₗ₀₀,
              Cᵇₖ₁   = Cᶜᵒᵉᶠᶠ.Cᵇₖ₁,
              Cᴴ²ᴼₖ₁ = Cᶜᵒᵉᶠᶠ.Cᴴ²ᴼₖ₁,
              Cᴮᵀ    = Cᶜᵒᵉᶠᶠ.Cᴮᵀ,
    )
    println("Borate concentration = ",Pᶜᵒᵉᶠᶠ.Cᴮᵀ)

    # Calculate pH, pCO2, and carbon species from Aᵀ and Cᵀ
    pH, CO₂ˢᵒˡ, HCO₃⁻, CO₃²⁻, _, pCO₂ᵒᶜᵉ, _  = Fᵖᶜᵒ²⁽ᴬᵀ⁺ᶜᵀ⁾(Aᵀ, Cᵀ, pH, Pᶜᵒᵉᶠᶠ)

    return CarbonSystem(pH, CO₂ˢᵒˡ, HCO₃⁻, CO₃²⁻, Cᵀ, Aᵀ, pCO₂ᵒᶜᵉ, pCO₂ᵃᵗᵐ)
end # end function

"""
    Fᶜᵀ⁽ᴬᵀ⁺ᵖᶜᵒ²⁾(Aᵀ, pCO₂ᵃᵗᵐ, pH, Pᶜᵒᵉᶠᶠ)

    Solve for DIC given total Alkalinity and atmosphere pCO₂
"""
@inline function Fᶜᵀ⁽ᴬᵀ⁺ᵖᶜᵒ²⁾(Aᵀ, pCO₂ᵃᵗᵐ, pH, Pᶜᵒᵉᶠᶠ)
    # Find the real roots of the polynomial using RootSolvers.jl 
    sol = find_zero(  x -> (
        x^3*(Aᵀ) +
        x^2*(
             Pᶜᵒᵉᶠᶠ.Cᵇₖ₁*
             Aᵀ-
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀*
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*
             pCO₂ᵃᵗᵐ-
             Pᶜᵒᵉᶠᶠ.Cᵇₖ₁*
             Pᶜᵒᵉᶠᶠ.Cᴮᵀ) + 
        x^1*(
            -Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀*
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*
             Pᶜᵒᵉᶠᶠ.Cᵇₖ₁*
             pCO₂ᵃᵗᵐ-
           2*Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀*
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂*
             pCO₂ᵃᵗᵐ
             ) +
        X^0*(
          -2*Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀*
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂*
             Pᶜᵒᵉᶠᶠ.Cᵇₖ₁*
             pCO₂ᵃᵗᵐ
             )        
        ),
        NewtonsMethodAD{Float64}(10^-(pH)),
        #SecantMethod{Float64}(10^-(pH+0.1), 10^-(pH-0.1)),
        CompactSolution());
    
    if sol.converged == true
        H = sol.root
        # Update pH
        pH = -log10(H)

        CO₂ˢᵒˡ = Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀*pCO₂ᵃᵗᵐ
        HCO₃⁻  = (Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*CO₂ˢᵒˡ)/H
        CO₃²⁻  = (Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂*CO₂ˢᵒˡ)/(H*H)

        return  pH, 
                CO₂ˢᵒˡ, 
                HCO₃⁻, 
                CO₃²⁻,
                CO₂ˢᵒˡ + HCO₃⁻ + CO₃²⁻, 
                CO₂ˢᵒˡ/Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀, 
                Aᵀ
    else
        error("CarbonSolverApprox did not converge")
        return nothing
    end
end # end function

"""
    Fᵖᶜᵒ²⁽ᴬᵀ⁺ᶜᵀ⁾(Aᵀ, Cᵀ, pH, Pᶜᵒᵉᶠᶠ)

    Solve for ocean pCO₂ given total Alkalinity and DIC
"""
@inline function Fᵖᶜᵒ²⁽ᴬᵀ⁺ᶜᵀ⁾(Aᵀ, Cᵀ, pH, Pᶜᵒᵉᶠᶠ)
    println("call to calculate pCO₂")
    # Find the real roots of the polynomial using RootSolvers.jl
    sol = find_zero(  x -> (
        x^3*(Aᵀ) +
        x^2*(
             Pᶜᵒᵉᶠᶠ.Cᵇₖ₁*
             Aᵀ
            +Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*
             Aᵀ
            -Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*
             Cᵀ+
            -Pᶜᵒᵉᶠᶠ.Cᵇₖ₁*
             Pᶜᵒᵉᶠᶠ.Cᴮᵀ
            ) + 
        x^1*(
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂*
             Aᵀ
          -2*Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂*
             Cᵀ
            +Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*
             Pᶜᵒᵉᶠᶠ.Cᵇₖ₁*
             Aᵀ 
            -Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*
             Pᶜᵒᵉᶠᶠ.Cᵇₖ₁*
             Cᵀ+
            -Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*
             Pᶜᵒᵉᶠᶠ.Cᵇₖ₁*
             Pᶜᵒᵉᶠᶠ.Cᴮᵀ
            ) +
        x^0*( 
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂*
             Pᶜᵒᵉᶠᶠ.Cᵇₖ₁*
             Aᵀ
          -2*Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂*
             Pᶜᵒᵉᶠᶠ.Cᵇₖ₁*
             Cᵀ
            -Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂*
             Pᶜᵒᵉᶠᶠ.Cᵇₖ₁*
             Pᶜᵒᵉᶠᶠ.Cᴮᵀ
            )                       
        ),
        #SecantMethod{Float64}(10^-(pH+0.1), 10^-(pH-0.1)),
        NewtonsMethodAD{Float64}(10^-(pH)),
        CompactSolution());

    if sol.converged == true
        H=sol.root

        # Update pH
        pH = -log10(H)

        CO₂ˢᵒˡ = (H*H*Cᵀ)/(H*H+H*Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁+Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂)
        HCO₃⁻  = (Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*CO₂ˢᵒˡ)/H
        CO₃²⁻  = (Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂*CO₂ˢᵒˡ)/(H*H)
        
        println("calculated pCO₂ = ", CO₂ˢᵒˡ/Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀ * 1e6)
        return pH, 
                CO₂ˢᵒˡ, 
                HCO₃⁻, 
                CO₃²⁻,
                CO₂ˢᵒˡ + HCO₃⁻ + CO₃²⁻, 
                CO₂ˢᵒˡ/Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀, 
                Aᵀ
    else
        error("CarbonSolverApprox did not converge")
        return nothing
    end
end # end function

end # module CarbonSolverApprox

# ------------------------------------------------------------

using .CarbonSystemSolvers.CarbonSolverApprox

include("carbon_chemistry_coefficients.jl")

## This should go in the testing suite, eventually.
Θᶜ      = 25.0
Sᴬ      = 35.0
Δpᵦₐᵣ   = 0.0
Cᵀ      = 2050e-6 # umol/kg to mol/kg
Aᵀ      = 2350e-6 # umol/kg to mol/kg
pCO₂ᵃᵗᵐ = 280e-6  # uatm to atm
pH      = 8.0
FT = Float64

Cᶜᵒᵉᶠᶠ = CarbonChemistryCoefficients(Θᶜ, Sᴬ, Δpᵦₐᵣ)
"""
#Check the values of calculated constants
println("Cᵈⁱᶜₖ₀ = ",       Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀     )
println("Cᵈⁱᶜₖ₁ᵣ₉₃ = ",    Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ᵣ₉₃  )
println("Cᵈⁱᶜₖ₂ᵣ₉₃ = ",    Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ᵣ₉₃  )
println("Cᵈⁱᶜₖ₁ₘ₉₅ = ",    Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ₘ₉₅  )
println("Cᵈⁱᶜₖ₂ₘ₉₅ = ",    Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ₘ₉₅  )
println("Cᵈⁱᶜₖ₁ₗ₀₀ = ",     Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ₗ₀₀  )
println("Cᵈⁱᶜₖ₂ₗ₀₀ = ",     Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ₗ₀₀  )
println("Cᵇₖ₁ = ",         Cᶜᵒᵉᶠᶠ.Cᵇₖ₁       )
println("Cᴴ²ᴼₖ₁ = ",       Cᶜᵒᵉᶠᶠ.Cᴴ²ᴼₖ₁     )
println("Cᴾᴼ⁴ₖ₁ = ",       Cᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁     )
println("Cᴾᴼ⁴ₖ₂ = ",       Cᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂     )
println("Cᴾᴼ⁴ₖ₃ = ",       Cᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₃     )
println("Cˢⁱᵗₖ₁ = ",       Cᶜᵒᵉᶠᶠ.Cˢⁱᵗₖ₁     )
println("Cᴴ²ˢₖ₁ = ",       Cᶜᵒᵉᶠᶠ.Cᴴ²ˢₖ₁     )
println("Cᴺᴴ⁴ₖ₁ = ",       Cᶜᵒᵉᶠᶠ.Cᴺᴴ⁴ₖ₁     )
println("Cᴴᶠᵦ₁ = ",        Cᶜᵒᵉᶠᶠ.Cᴴᶠᵦ₁     )
println("Cᴴᶠₖ₁ = ",        Cᶜᵒᵉᶠᶠ.Cᴴᶠₖ₁      )
println("Cᴴˢᴼ⁴ₖ₁ = ",      Cᶜᵒᵉᶠᶠ.Cᴴˢᴼ⁴ₖ₁    )
println("Cᶜᵃˡᶜⁱᵗᵉₛₚ = ",   Cᶜᵒᵉᶠᶠ.Cᶜᵃˡᶜⁱᵗᵉₛₚ  )
println("Cᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ = ", Cᶜᵒᵉᶠᶠ.Cᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ)
println("Cᴮᵀ = " ,        Cᶜᵒᵉᶠᶠ.Cᴮᵀ         )
println("Cᶠᵀ = " ,        Cᶜᵒᵉᶠᶠ.Cᶠᵀ         )
println("Cᶜᵃ = ",         Cᶜᵒᵉᶠᶠ.Cᶜᵃ         )
println("Cˢᴼ⁴ = ",        Cᶜᵒᵉᶠᶠ.Cˢᴼ⁴        )
"""

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

(; pH, CO₂ˢᵒˡ, HCO₃⁻, CO₃²⁻, Cᵀ, Aᵀ, pCO₂ᵒᶜᵉ, pCO₂ᵃᵗᵐ) = 
CarbonSystemApprox(
        Θᶜ, Sᴬ, Δpᵦₐᵣ, Cᵀ, Aᵀ, pH, pCO₂ᵃᵗᵐ,
        )
        
println("Cᵀ = ", Cᵀ * 1e6  )
println("Aᵀ = ", Aᵀ * 1e6  )
println("pH = " , pH       )
println("pCO₂ᵃᵗᵐ = ", pCO₂ᵃᵗᵐ * 1e6)
println("pCO₂ᵒᶜᵉ = ", pCO₂ᵒᶜᵉ * 1e6)
end