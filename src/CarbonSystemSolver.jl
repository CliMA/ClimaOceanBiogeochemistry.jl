module CarbonSystemSolver

module CarbonSystemApprox
export CarbonSolverApprox

using DissociationConstants
using RootSolvers

struct CarbonSolverApprox
    "Temperature in degrees Celsius"
    Θᶜ    :: Float64
    "Absolute Salinity in g/kg"
    Sᴬ    :: Float64
    "Applied Pressure in bars"
    Δpᵦₐᵣ :: Float64
    "DIC concentration in mol C/kg"
    Cᵀ     :: Float64
    "Alkalinity in mol Eq/kg"
    Aᵀ    :: Float64
    "Atmospheric pCO₂ in atm"
    pCO₂ᵃᵗᵐ :: Float64
    """
    CarbonSolverApprox(
        FT            = Float64,
        Θ       :: FT = 25.0,
        Sᴬ      :: FT = 35.0,
        Δpᵦₐᵣ   :: FT = 0.0,
        Cᵀ      :: FT = 2050.0e-6,
        Aᵀ      :: FT = 2350.0e-6,
        pCO₂ᵃᵗᵐ :: FT = 280.0e-6,
        )

TBW
"""
function CarbonSolverApprox(
        Θᶜ      :: Float64 = 25.0,
        Sᴬ      :: Float64 = 35.0,
        Δpᵦₐᵣ   :: Float64 = 0.0,
        Cᵀ      :: Float64 = 2050.0e-6,
        Aᵀ      :: Float64 = 2350.0e-6,
        pH      :: Float64 = 8.0,
        pCO₂ᵃᵗᵐ :: Float64 = 280.0e-6,
        )

        # DissociationCoefficients are pretty much all in mol/kg, hence the 1e-6 factors for Cᵀ and Aᵀ
        Cᶜᵒⁿˢᵗ = DissociationCoefficients(Θᶜ, Sᴬ, Δpᵦₐᵣ)

        # Some logic here about choosing coefficient options, particularly Cᵈⁱᶜ 
        Pᶜᵒⁿˢᵗ = (Cᵈⁱᶜₖ₀ = Cᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₀,
                  Cᵈⁱᶜₖ₁ = Cᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁ₗ₀₀,
                  Cᵈⁱᶜₖ₂ = Cᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₂ₗ₀₀,
                  Cᵇₖ₁   = Cᶜᵒⁿˢᵗ.Cᵇₖ₁,
                  Cᴴ²ᴼₖ₁ = Cᶜᵒⁿˢᵗ.Cᴴ²ᴼₖ₁,
                  Cᴮᵀ    = Cᶜᵒⁿˢᵗ.Cᴮᵀ,
        )
        return Fᵖᶜᵒ²⁽ᴬᵀ⁺ᶜᵀ⁾(Aᵀ,Cᵀ,pH,Pᶜᵒⁿˢᵗ)
    end # end function
end # end struct

@inline """
    Fᶜᵀ⁽ᴬᵀ⁺ᵖᶜᵒ²⁾(Aᵀ,pCO₂ᵃᵗᵐ,pH,Pᶜᵒⁿˢᵗ)

    Solve for DIC given total Alkalinity and atmosphere pCO₂
"""
function Fᶜᵀ⁽ᴬᵀ⁺ᵖᶜᵒ²⁾(Aᵀ,pCO₂ᵃᵗᵐ,pH,Pᶜᵒⁿˢᵗ)
    # Find the real roots of the polynomial using RootSolvers.jl 
    sol = find_zero(  x -> (
        x^3*(Aᵀ) +
        x^2*(
             Pᶜᵒⁿˢᵗ.Cᵇₖ₁*
             Aᵀ-
             Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₀*
             Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁*
             pCO₂ᵃᵗᵐ-
             Pᶜᵒⁿˢᵗ.Cᵇₖ₁*
             Pᶜᵒⁿˢᵗ.Cᴮᵀ) + 
        x^1*(
            -Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₀*
             Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁*
             Pᶜᵒⁿˢᵗ.Cᵇₖ₁*
             pCO₂ᵃᵗᵐ-
           2*Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₀*
             Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁*
             Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₂*
             pCO₂ᵃᵗᵐ
             ) +
        X^0*(
          -2*Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₀*
             Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁*
             Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₂*
             Pᶜᵒⁿˢᵗ.Cᵇₖ₁*
             pCO₂ᵃᵗᵐ
             )        
        ),
        #NewtonsMethodAD{Float64}(10^-pH),
        SecantMethod{Float64}(10^-(pH+1.0), 10^-(pH-1.0)),
        CompactSolution());
    
    if sol.converged == true
        H = sol.root
        # Update pH
        pH = -log10(H)

        CO₂ˢᵒˡ = Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₀*pCO₂ᵃᵗᵐ
        HCO₃⁻  = (Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁*CO₂ˢᵒˡ)/H
        CO₃²⁻  = (Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁*Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₂*CO₂ˢᵒˡ)/(H*H)

        return CO₂ˢᵒˡ + HCO₃⁻ + CO₃²⁻, pH
    else
        error("CarbonSolverApprox did not converge")
        return NaN
    end
end

@inline """
    Fᵖᶜᵒ²⁽ᴬᵀ⁺ᶜᵀ⁾(Aᵀ,Cᵀ,pH,Pᶜᵒⁿˢᵗ)

    Solve for pCO₂ given total Alkalinity and DIC
"""
function Fᵖᶜᵒ²⁽ᴬᵀ⁺ᶜᵀ⁾(Aᵀ,Cᵀ,pH,Pᶜᵒⁿˢᵗ)
    # Find the real roots of the polynomial using RootSolvers.jl
    sol = find_zero(  x -> (
        x^3*(Aᵀ) +
        x^2*(
             Pᶜᵒⁿˢᵗ.Cᵇₖ₁*
             Aᵀ
            +Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁*
             Aᵀ
            -Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁*
             Cᵀ+
            -Pᶜᵒⁿˢᵗ.Cᵇₖ₁*
             Pᶜᵒⁿˢᵗ.Cᴮᵀ
            ) + 
        x^1*(
             Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁*
             Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₂*
             Aᵀ
          -2*Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁*
             Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₂*
             Cᵀ
            +Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁*
             Pᶜᵒⁿˢᵗ.Cᵇₖ₁*
             Aᵀ 
            -Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁*
             Pᶜᵒⁿˢᵗ.Cᵇₖ₁*
             Cᵀ+
            -Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁*
             Pᶜᵒⁿˢᵗ.Cᵇₖ₁*
             Pᶜᵒⁿˢᵗ.Cᴮᵀ
            ) +
        x^0*( 
             Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁*
             Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₂*
             Pᶜᵒⁿˢᵗ.Cᵇₖ₁*
             Aᵀ
          -2*Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁*
             Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₂*
             Pᶜᵒⁿˢᵗ.Cᵇₖ₁*
             Cᵀ
            -Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁*
             Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₂*
             Pᶜᵒⁿˢᵗ.Cᵇₖ₁*
             Pᶜᵒⁿˢᵗ.Cᴮᵀ
            )                       
        ),
        NewtonsMethodAD{Float64}(10^-pH),
        #SecantMethod{Float64}(10^-(pH+1.0), 10^-(pH-1.0)),
        CompactSolution());

    if sol.converged == true
        H=sol.root

        # Update pH
        pH = -log10(H)

        CO₂ˢᵒˡ = (H*H*Cᵀ)/(H*H+H*Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁+Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁*Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₂)
        HCO₃⁻  = (Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁*CO₂ˢᵒˡ)/H
        CO₃²⁻  = (Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁*Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₂*CO₂ˢᵒˡ)/(H*H)
        return CO₂ˢᵒˡ/Pᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₀, pH
    else
        error("CarbonSolverApprox did not converge")
        return NaN
    end
end

end # module CarbonSolverApprox
# ------------------------------------------------------------

using .CarbonSystemApprox
using DissociationConstants

## This should go in the testing suite, eventually.
Θᶜ      = 25.0
Sᴬ      = 35.0
Δpᵦₐᵣ   = 0.0
Cᵀ      = 2050.0*1e-6 # umol/kg to mol/kg
Aᵀ      = 2350.0*1e-6 # umol/kg to mol/kg
pCO₂ᵃᵗᵐ = 280.0*1e-6  # uatm to atm
pH      = 8.0

Cᶜᵒⁿˢᵗ = DissociationCoefficients(Θᶜ, Sᴬ, Δpᵦₐᵣ);
"""
#Check the values of calculated constants
println("Cᵈⁱᶜₖ₀ = ",       Cᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₀     )
println("Cᵈⁱᶜₖ₁ᵣ₉₃ = ",    Cᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁ᵣ₉₃  )
println("Cᵈⁱᶜₖ₂ᵣ₉₃ = ",    Cᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₂ᵣ₉₃  )
println("Cᵈⁱᶜₖ₁ₘ₉₅ = ",    Cᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁ₘ₉₅  )
println("Cᵈⁱᶜₖ₂ₘ₉₅ = ",    Cᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₂ₘ₉₅  )
println("Cᵈⁱᶜₖ₁ₗ₀₀ = ",     Cᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁ₗ₀₀  )
println("Cᵈⁱᶜₖ₂ₗ₀₀ = ",     Cᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₂ₗ₀₀  )
println("Cᵇₖ₁ = ",         Cᶜᵒⁿˢᵗ.Cᵇₖ₁       )
println("Cᴴ²ᴼₖ₁ = ",       Cᶜᵒⁿˢᵗ.Cᴴ²ᴼₖ₁     )
println("Cᴾᴼ⁴ₖ₁ = ",       Cᶜᵒⁿˢᵗ.Cᴾᴼ⁴ₖ₁     )
println("Cᴾᴼ⁴ₖ₂ = ",       Cᶜᵒⁿˢᵗ.Cᴾᴼ⁴ₖ₂     )
println("Cᴾᴼ⁴ₖ₃ = ",       Cᶜᵒⁿˢᵗ.Cᴾᴼ⁴ₖ₃     )
println("Cˢⁱᵗₖ₁ = ",       Cᶜᵒⁿˢᵗ.Cˢⁱᵗₖ₁     )
println("Cᴴ²ˢₖ₁ = ",       Cᶜᵒⁿˢᵗ.Cᴴ²ˢₖ₁     )
println("Cᴺᴴ⁴ₖ₁ = ",       Cᶜᵒⁿˢᵗ.Cᴺᴴ⁴ₖ₁     )
println("Cᴴᶠᵦ₁ = ",        Cᶜᵒⁿˢᵗ.Cᴴᶠᵦ₁     )
println("Cᴴᶠₖ₁ = ",        Cᶜᵒⁿˢᵗ.Cᴴᶠₖ₁      )
println("Cᴴˢᴼ⁴ₖ₁ = ",      Cᶜᵒⁿˢᵗ.Cᴴˢᴼ⁴ₖ₁    )
println("Cᶜᵃˡᶜⁱᵗᵉₛₚ = ",   Cᶜᵒⁿˢᵗ.Cᶜᵃˡᶜⁱᵗᵉₛₚ  )
println("Cᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ = ", Cᶜᵒⁿˢᵗ.Cᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ)
println("Cᴮᵀ = " ,        Cᶜᵒⁿˢᵗ.Cᴮᵀ         )
println("Cᶠᵀ = " ,        Cᶜᵒⁿˢᵗ.Cᶠᵀ         )
println("Cᶜᵃ = ",         Cᶜᵒⁿˢᵗ.Cᶜᵃ         )
println("Cˢᴼ⁴ = ",        Cᶜᵒⁿˢᵗ.Cˢᴼ⁴        )
"""
@assert round(log(Cᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₀), digits = 4)     == -3.5617 # Handbook (2007)
@assert round(log(Cᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁ᵣ₉₃), digits = 4)  == -13.4847 # Handbook (1994)
@assert round(log(Cᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₂ᵣ₉₃), digits = 4)  == -20.5504 # Handbook (1994)
#@assert Cᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁ₘ₉₅  ==
#@assert Cᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₂ₘ₉₅  ==
@assert round(log10(Cᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₁ₗ₀₀), digits = 4)   == -5.8472 # Handbook (2007)
@assert round(log10(Cᶜᵒⁿˢᵗ.Cᵈⁱᶜₖ₂ₗ₀₀), digits = 4)   == -8.9660 # Handbook (2007)
@assert round(log(Cᶜᵒⁿˢᵗ.Cᵇₖ₁), digits = 4)       == -19.7964 # Handbook (2007)
@assert round(log(Cᶜᵒⁿˢᵗ.Cᴴ²ᴼₖ₁*Cᶜᵒⁿˢᵗ.CH⁺ₛoverH⁺ₜ)-0.015, digits = 3)     == -30.434 # Handbook (2007)
@assert round(log(Cᶜᵒⁿˢᵗ.Cᴾᴼ⁴ₖ₁*Cᶜᵒⁿˢᵗ.CH⁺ₛoverH⁺ₜ)-0.015, digits = 2)     == -3.71 # Handbook (2007)
@assert round(log(Cᶜᵒⁿˢᵗ.Cᴾᴼ⁴ₖ₂*Cᶜᵒⁿˢᵗ.CH⁺ₛoverH⁺ₜ)-0.015, digits = 3)     == -13.727 # Handbook (2007)
@assert round(log(Cᶜᵒⁿˢᵗ.Cᴾᴼ⁴ₖ₃*Cᶜᵒⁿˢᵗ.CH⁺ₛoverH⁺ₜ)-0.015, digits = 2)     == -20.24 # Handbook (2007)
@assert round(log(Cᶜᵒⁿˢᵗ.Cˢⁱᵗₖ₁*Cᶜᵒⁿˢᵗ.CH⁺ₛoverH⁺ₜ)-0.015, digits = 2)     == -21.61 # Handbook (2007)
@assert round(-log10(Cᶜᵒⁿˢᵗ.Cᴴ²ˢₖ₁*Cᶜᵒⁿˢᵗ.CH⁺ₛoverH⁺ₜ), digits = 2)     == 6.51 # Lewis and Wallace (1998)
@assert round(-log10(Cᶜᵒⁿˢᵗ.Cᴺᴴ⁴ₖ₁*Cᶜᵒⁿˢᵗ.CH⁺ₛoverH⁺ₜ), digits = 2)     == 9.26 # Lewis and Wallace (1998)
#@assert Cᶜᵒⁿˢᵗ.Cᴴᶠᵦ₁      ==
@assert round(log(Cᶜᵒⁿˢᵗ.Cᴴᶠₖ₁), digits = 2)       == -6.09 # Handbook (2007)
@assert round(log(Cᶜᵒⁿˢᵗ.Cᴴˢᴼ⁴ₖ₁), digits = 2)     == -2.30 # Handbook (2007)
#@assert Cᶜᵒⁿˢᵗ.Cᶜᵃˡᶜⁱᵗᵉₛₚ  ==
#@assert Cᶜᵒⁿˢᵗ.Cᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ==
#@assert Cᶜᵒⁿˢᵗ.Cᴮᵀ        ==
#@assert Cᶜᵒⁿˢᵗ.Cᶠᵀ        ==
#@assert Cᶜᵒⁿˢᵗ.Cᶜᵃ        ==
#@assert Cᶜᵒⁿˢᵗ.Cˢᴼ⁴       ==

#pCO₂, pH = CarbonSolverApprox(
#        Θᶜ, Sᴬ, Δpᵦₐᵣ, Cᵀ, Aᵀ, pH, pCO₂ᵃᵗᵐ,
#        )
#println("Cᵀ = ", Cᵀ*1e6    )
#println("Aᵀ = ", Aᵀ*1e6    )
#println("pH = "  , pH      )
#println("pCO₂ᵃᵗᵐ = ", pCO₂ᵃᵗᵐ*1e6)
#println("pCO₂ = ", pCO₂*1e6)

end # module CarbonSolver