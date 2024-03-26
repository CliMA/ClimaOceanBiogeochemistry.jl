module CarbonSystemSolvers
export CarbonSystem

include("carbon_chemistry_coefficients.jl")

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
    FCᵀCO₂ˢᵒˡ(Cᵀ, pH, Pᶜᵒᵉᶠᶠ)

Calculate the dissolved and hydrated CO₂ concentration in seawater
given the total carbon concentration Cᵀ, pH, and the carbon chemistry coefficients.
"""
@inline function FCᵀCO₂ˢᵒˡ(Cᵀ, pH, Pᶜᵒᵉᶠᶠ)
    # Calculate H⁺ from pH
    H⁺ = 10^-pH

    return ( Cᵀ * H⁺^2 )/(
             (H⁺^2) + 
             (H⁺^1 * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁) + 
             (H⁺^0 * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂)
            )
end

"""
    FCᵀHCO₃⁻(Cᵀ, pH, Pᶜᵒᵉᶠᶠ)

Calculate the bicarbonate ion concentration in seawater
given the total carbon concentration Cᵀ, pH, and the carbon chemistry coefficients.
"""
@inline function FCᵀHCO₃⁻(Cᵀ, pH, Pᶜᵒᵉᶠᶠ)
    # Calculate H⁺ from pH
    H⁺ = 10^-pH

    return ( Cᵀ * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ * H⁺ )/(
             (H⁺^2) + 
             (H⁺^1 * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁) + 
             (H⁺^0 * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂)
            )
end

"""
    FCᵀCO₃²⁻(Cᵀ, pH, Pᶜᵒᵉᶠᶠ)

Calculate the carbonate concentration in seawater
given the total carbon concentration Cᵀ, pH, and the carbon chemistry coefficients.
"""
@inline function FCᵀCO₃²⁻(Cᵀ, pH, Pᶜᵒᵉᶠᶠ)
    # Calculate H⁺ from pH
    H⁺ = 10^-pH

    return ( Cᵀ * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ ) / (
             (H⁺^2) + 
             (H⁺^1 * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁) + 
             (H⁺^0 * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂)
            )
end

"""
    FpCO₂CO₂ˢᵒˡ(pCO₂, pH, Pᶜᵒᵉᶠᶠ)

Calculate the dissolved and hydrated CO₂ concentration in seawater
given the pCO₂, pH, and the carbon chemistry coefficients.
"""
@inline function FpCO₂CO₂ˢᵒˡ(pCO₂, Pᶜᵒᵉᶠᶠ)
    # Perhaps take account of fugacity here?
    return Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀ * pCO₂
end

"""
    FpCO₂HCO₃⁻(pCO₂, pH, Pᶜᵒᵉᶠᶠ)

Calculate the bicarbonate ion concentration in seawater
given the pCO₂, pH, and the carbon chemistry coefficients.
"""
@inline function FpCO₂HCO₃⁻(pCO₂, pH, Pᶜᵒᵉᶠᶠ)
    # Calculate H⁺ from pH
    H⁺ = 10^-pH

    return (Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀ * pCO₂)/H⁺
end

"""
    FpCO₂CO₃²⁻(pCO₂, pH, Pᶜᵒᵉᶠᶠ)

Calculate the carbonate concentration in seawater
given the pCO₂, pH, and the carbon chemistry coefficients.
"""
@inline function FpCO₂CO₃²⁻(pCO₂, pH, Pᶜᵒᵉᶠᶠ)
    # Calculate H⁺ from pH
    H⁺ = 10^-pH

    return (Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀ * pCO₂)/(H⁺^2)
end

# ----------------------------------------------------------------------------------
module UniversalRobustCarbonSolver
export UniversalRobustCarbonSystem,
        CarbonSystem

using ..CarbonSystemSolvers: CarbonSystem, CarbonChemistryCoefficients, FCᵀCO₂ˢᵒˡ, FCᵀCO₃²⁻, FCᵀHCO₃⁻
#include("carbon_chemistry_coefficients.jl")

"""
    UniversalRobustCarbonSystem(
            Θ       :: FT = 25.0,
            Sᴬ      :: FT = 35.0,
            Δpᵦₐᵣ   :: FT = 0.0,
            Cᵀ      :: FT = 2050.0e-6,
            Aᵀ      :: FT = 2350.0e-6,
            Pᵀ      :: FT = 1.0e-6,
            Siᵀ     :: FT = 15.0e-6,
            pH      :: FT = 8.0,
            pCO₂ᵃᵗᵐ :: FT = 280.0e-6,
            )

Uses the Munhoven (2013) SolveSAPHE package to solve the distribution of carbon species
"""
@inline function UniversalRobustCarbonSystem(
        Θᶜ      :: FT = 25.0,
        Sᴬ      :: FT = 35.0,
        Δpᵦₐᵣ   :: FT = 0.0,
        Cᵀ      :: FT = 2050.0e-6,
        Aᵀ      :: FT = 2350.0e-6,
        Pᵀ      :: FT = 1.0e-6,
        Siᵀ     :: FT = 15.0e-6,
        pH      :: FT = 8.0,
        pCO₂ᵃᵗᵐ :: FT = 280.0e-6) where {FT}

    # CarbonChemistryCoefficients are pretty much all in mol/kg, hence the 1e-6 factors for Cᵀ and Aᵀ
    Cᶜᵒᵉᶠᶠ = CarbonChemistryCoefficients(Θᶜ, Sᴬ, Δpᵦₐᵣ)
    
    # Some logic here about choosing coefficient options, particularly Cᵈⁱᶜ 
    Pᶜᵒᵉᶠᶠ = (Cᵈⁱᶜₖ₀ = Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀,
              Cᵈⁱᶜₖ₁ = Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ₗ₀₀,
              Cᵈⁱᶜₖ₂ = Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ₗ₀₀,
              Cᵇₖ₁   = Cᶜᵒᵉᶠᶠ.Cᵇₖ₁,
              Cᴾᴼ⁴ₖ₁ = Cᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁,
              Cᴾᴼ⁴ₖ₂ = Cᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂,
              Cᴾᴼ⁴ₖ₃ = Cᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₃,
              Cˢⁱᵗₖ₁ = Cᶜᵒᵉᶠᶠ.Cˢⁱᵗₖ₁,
              Cᴴˢᴼ⁴ₖ₁= Cᶜᵒᵉᶠᶠ.Cᴴˢᴼ⁴ₖ₁,
              Cᴴ²ˢₖ₁ = Cᶜᵒᵉᶠᶠ.Cᴴ²ˢₖ₁,
              Cᴴᶠₖ₁  = Cᶜᵒᵉᶠᶠ.Cᴴᶠₖ₁,
              Cᴺᴴ⁴ₖ₁ = Cᶜᵒᵉᶠᶠ.Cᴺᴴ⁴ₖ₁,
              Cᴴ²ᴼₖ₁ = Cᶜᵒᵉᶠᶠ.Cᴴ²ᴼₖ₁,
              Cᴮᵀ    = Cᶜᵒᵉᶠᶠ.Cᴮᵀ,
              Cᶠᵀ    = Cᶜᵒᵉᶠᶠ.Cᶠᵀ,
              Cˢᴼ⁴   = Cᶜᵒᵉᶠᶠ.Cˢᴼ⁴,
              H⁺ₜoverH⁺₃ = Cᶜᵒᵉᶠᶠ.H⁺ₜoverH⁺₃, # pH scale conversion factor
    )

    # Calculate pH from Aᵀ and Cᵀ and then calculate the rest of the carbon system
    pH     = Fᵖᴴᵤₙᵢᵣₒ(Aᵀ, Cᵀ, Pᵀ, Siᵀ, pH, Pᶜᵒᵉᶠᶠ) 
    CO₂ˢᵒˡ = FCᵀCO₂ˢᵒˡ(Cᵀ, pH, Pᶜᵒᵉᶠᶠ)
    HCO₃⁻  = FCᵀHCO₃⁻(Cᵀ, pH, Pᶜᵒᵉᶠᶠ)
    CO₃²⁻  = FCᵀCO₃²⁻(Cᵀ, pH, Pᶜᵒᵉᶠᶠ)
    pCO₂ᵒᶜᵉ= CO₂ˢᵒˡ / Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀ # correct for fugacity of CO₂ in seawater?

    return CarbonSystem(pH, CO₂ˢᵒˡ, HCO₃⁻, CO₃²⁻, Cᵀ, Aᵀ, pCO₂ᵒᶜᵉ, pCO₂ᵃᵗᵐ)
end # end function

"""
    Fᵖᴴᵤₙᵢᵣₒ(Aᵀ, Cᵀ, Pᵀ, Siᵀ, pH, Pᶜᵒᵉᶠᶠ)

Calculate the pH of seawater given the total alkalinity Aᵀ, total carbon Cᵀ,
total phosphate Pᵀ, total silicate Siᵀ, and the carbon chemistry coefficients.
Uses the SolveSAPHE package (Munhoven et al., 2013), a universal, robust, pH 
solver that converges from any given initial value.
"""
@inline function Fᵖᴴᵤₙᵢᵣₒ(Aᵀ, Cᵀ, Pᵀ, Siᵀ, pH, Pᶜᵒᵉᶠᶠ, NH₄ᵀ=0, H₂Sᵀ=0, Δₕ₊=1e-8, eᴴ⁺ᵗʰʳᵉˢʰ=1, Iᴴ⁺ₘₐₓ=100) 
   
    Iᴴ⁺                = 0
    Aᵀᵃᵇˢₘᵢₙ           = Inf
    H⁺ᶠᵃᶜᵗᵒʳ           = 1

    if pH == 8
        # Get a better initial H+ guess
        H⁺ᵢₙᵢ = FH⁺ᵢₙᵢ(Aᵀ, Cᵀ, Pᶜᵒᵉᶠᶠ)
    else
        # Calculate H⁺ from pH
        H⁺ᵢₙᵢ = 10^-pH
    end

    # Calculate initial bounds of H+ concentration
    Aᵀₗₒ, Aᵀₕᵢ = FboundsAᵀₙₕ₂ₒ( 
                    Cᵀ, Pᵀ, Siᵀ, NH₄ᵀ, H₂Sᵀ, Pᶜᵒᵉᶠᶠ
                  )

    Δₗₒ = (Aᵀ - Aᵀₗₒ)^2 + 4 * Pᶜᵒᵉᶠᶠ.Cᴴ²ᴼₖ₁/Pᶜᵒᵉᶠᶠ.H⁺ₜoverH⁺₃

    if Aᵀ ≥ Aᵀₗₒ
        H⁺ₘᵢₙ = 2 * Pᶜᵒᵉᶠᶠ.Cᴴ²ᴼₖ₁/( Aᵀ - Aᵀₗₒ + sqrt(Δₗₒ) )
    else
        H⁺ₘᵢₙ = Pᶜᵒᵉᶠᶠ.H⁺ₜoverH⁺₃ *( -(Aᵀ - Aᵀₗₒ) + sqrt(Δₗₒ) )/2
    end

    Δₕᵢ = (Aᵀ - Aᵀₕᵢ)^2 + 4 * Pᶜᵒᵉᶠᶠ.Cᴴ²ᴼₖ₁/Pᶜᵒᵉᶠᶠ.H⁺ₜoverH⁺₃

    if Aᵀ ≤ Aᵀₕᵢ
        H⁺ₘₐₓ = Pᶜᵒᵉᶠᶠ.H⁺ₜoverH⁺₃ *( -(Aᵀ - Aᵀₕᵢ) + sqrt(Δₕᵢ) )/2
    else
        H⁺ₘₐₓ = 2 * Pᶜᵒᵉᶠᶠ.Cᴴ²ᴼₖ₁/( Aᵀ - Aᵀₕᵢ + sqrt(Δₕᵢ) )
    end

    # Initial guess for H⁺
    H⁺ = max(min(H⁺ₘₐₓ, H⁺ᵢₙᵢ), H⁺ₘᵢₙ)             
    #H⁺ = sqrt(H⁺ₘₐₓ * H⁺ₘᵢₙ) # Safer(?) than the above line

    while abs(H⁺ᶠᵃᶜᵗᵒʳ) > Δₕ₊
    # Stop iterations once |\delta{[H]}/[H]| < rdel
    # <=> |(H⁺ - H⁺ₚᵣₑ)/H⁺ₚᵣₑ| = |EXP(-Aᵀᵣₐₜ/(∂Aᵀᵣₐₜ∂H⁺*H⁺ₚᵣₑ)) -1| < rdel
    # |EXP(-Aᵀᵣₐₜ/(∂Aᵀᵣₐₜ∂H⁺*H⁺ₚᵣₑ)) -1| ~ |Aᵀᵣₐₜ/(∂Aᵀᵣₐₜ∂H⁺*H⁺ₚᵣₑ)|
    # Alternatively:
    # |\Delta pH| = |Aᵀᵣₐₜ/(∂Aᵀᵣₐₜ∂H⁺*H⁺ₚᵣₑ*LOG(10))|
    #             ~ 1/LOG(10) * |\Delta [H]|/[H]
    #             < 1/LOG(10) * rdel
    # Hence |Aᵀᵣₐₜ/(∂Aᵀᵣₐₜ∂H⁺*H⁺)| < rdel
    # rdel <-- Δₕ₊
     
        if Iᴴ⁺ ≥ Iᴴ⁺ₘₐₓ
            H⁺ = nothing
            break
        end
         
        # Increase the iteration counter
        Iᴴ⁺ += 1

        # remember for next iteration current H⁺ concentration
        H⁺ₚᵣₑ = H⁺

        Aᵀᵣₐₜ, ∂Aᵀᵣₐₜ∂H⁺ = FAᵀ(
                              Cᵀ, Aᵀ, Pᵀ, Siᵀ, NH₄ᵀ, H₂Sᵀ, H⁺, Pᶜᵒᵉᶠᶠ
                             )

        # Adapt bracketing interval
        if Aᵀᵣₐₜ > 0
           H⁺ₘᵢₙ = H⁺ₚᵣₑ
        elseif Aᵀᵣₐₜ < 0
           H⁺ₘₐₓ = H⁺ₚᵣₑ
        else
        # H⁺ is the root; unlikely but, one never knows
           break
        end

        if abs(Aᵀᵣₐₜ) ≥ Aᵀᵃᵇˢₘᵢₙ/2
        # if the function evaluation at the current point is
        # not decreasing faster than with a bisection step (at least linearly)
        # in absolute value take one bisection step on [ph_min, ph_max]
        # ph_new = (ph_min + ph_max)/2d0
        # In terms of [H]_new:
        # [H]_new = 10**(-ph_new)
        #         = 10**(-(ph_min + ph_max)/2d0)
        #         = SQRT(10**(-(ph_min + phmax)))
        #         = SQRT(H⁺ₘₐₓ * H⁺ₘᵢₙ)
     
            H⁺        = sqrt(H⁺ₘₐₓ * H⁺ₘᵢₙ)
            H⁺ᶠᵃᶜᵗᵒʳ  = ( H⁺ - H⁺ₚᵣₑ ) / H⁺ₚᵣₑ 
        else
        # dAᵀᵣₐₜ/dpH = dAᵀᵣₐₜ/d[H] * d[H]/dpH
        #           = -∂Aᵀᵣₐₜ∂H⁺ * LOG(10) * [H]
        # \Delta pH = -Aᵀᵣₐₜ/(∂Aᵀᵣₐₜ∂H⁺*d[H]/dpH) = Aᵀᵣₐₜ/(∂Aᵀᵣₐₜ∂H⁺*[H]*LOG(10))
        # pH_new = pH_old + \deltapH
        # [H]_new = 10**(-pH_new)
        #         = 10**(-pH_old - \Delta pH)
        #         = [H]_old * 10**(-Aᵀᵣₐₜ/(∂Aᵀᵣₐₜ∂H⁺*[H]_old*LOG(10)))
        #         = [H]_old * EXP(-LOG(10)*Aᵀᵣₐₜ/(∂Aᵀᵣₐₜ∂H⁺*[H]_old*LOG(10)))
        #         = [H]_old * EXP(-Aᵀᵣₐₜ/(∂Aᵀᵣₐₜ∂H⁺*[H]_old))

            H⁺ᶠᵃᶜᵗᵒʳ = -Aᵀᵣₐₜ / ( ∂Aᵀᵣₐₜ∂H⁺ * H⁺ₚᵣₑ )

            if abs(H⁺ᶠᵃᶜᵗᵒʳ) > eᴴ⁺ᵗʰʳᵉˢʰ
               H⁺ = H⁺ₚᵣₑ * exp( H⁺ᶠᵃᶜᵗᵒʳ )
            else
               H⁺ = H⁺ₚᵣₑ + (H⁺ᶠᵃᶜᵗᵒʳ * H⁺ₚᵣₑ)
            end

            if H⁺ < H⁺ₘᵢₙ
            # if [H]_new < [H]_min
            # i.e., if ph_new > ph_max then
            # take one bisection step on [ph_prev, ph_max]
            # ph_new = (ph_prev + ph_max)/2d0
            # In terms of [H]_new:
            # [H]_new = 10**(-ph_new)
            #         = 10**(-(ph_prev + ph_max)/2d0)
            #         = SQRT(10**(-(ph_prev + phmax)))
            #         = SQRT([H]_old*10**(-ph_max))
            #         = SQRT([H]_old * H⁺ₘᵢₙ)
               H⁺        = sqrt( H⁺ₚᵣₑ * H⁺ₘᵢₙ )
               H⁺ᶠᵃᶜᵗᵒʳ  = ( H⁺ - H⁺ₚᵣₑ ) / H⁺ₚᵣₑ 
            end

            if H⁺ > H⁺ₘₐₓ 
            # if [H]_new > [H]_max
            # i.e., if ph_new < ph_min, then
            # take one bisection step on [ph_min, ph_prev]
            # ph_new = (ph_prev + ph_min)/2d0
            # In terms of [H]_new:
            # [H]_new = 10**(-ph_new)
            #         = 10**(-(ph_prev + ph_min)/2d0)
            #         = SQRT(10**(-(ph_prev + ph_min)))
            #         = SQRT([H]_old*10**(-ph_min))
            #         = SQRT([H]_old * zhmax)
               H⁺       = sqrt( H⁺ₚᵣₑ * H⁺ₘₐₓ )
               H⁺ᶠᵃᶜᵗᵒʳ = ( H⁺ - H⁺ₚᵣₑ ) / H⁺ₚᵣₑ
            end
        end

        Aᵀᵃᵇˢₘᵢₙ = min( abs(Aᵀᵣₐₜ), Aᵀᵃᵇˢₘᵢₙ)
    end # end while loop

    if H⁺ > 0
        Aᵀᵣₐₜ, ∂Aᵀᵣₐₜ∂H⁺ = FAᵀ(
            Cᵀ, Aᵀ, Pᵀ, Siᵀ, NH₄ᵀ, H₂Sᵀ, H⁺, Pᶜᵒᵉᶠᶠ
        )
    else
        ∂Aᵀᵣₐₜ∂H⁺ = nothing
    end
    return -log10(H⁺)
end

"""
    FH⁺ᵢₙᵢ(Aᵀ, Cᵀ, H⁺, Pᶜᵒᵉᶠᶠ)

Calculates the root for the 2nd order approximation of the
 Cᵀ-Bᵀ-Aᶜ equation for [H+] (reformulated as a cubic polynomial)
 around the local minimum, if it exists.

 Returns * 1e-03 if Aᶜ <= 0
         * 1e-10 if Aᶜ >= 2*Cᵀ + Bᵀ
         * 1e-07 if 0 < Aᶜ < 2*Cᵀ + Bᵀ
            and the 2nd order approximation does not have a solution

"""
@inline function FH⁺ᵢₙᵢ(Aᶜ, Cᵀ, Pᶜᵒᵉᶠᶠ)

    if Aᶜ <= 0
        return 1e-3
    elseif Aᶜ >= (2*Cᵀ + Pᶜᵒᵉᶠᶠ.Cᴮᵀ)
        return 1e-10
    else
        Rᶜᴬ = Cᵀ/Aᶜ
        Rᴮᴬ = Pᶜᵒᵉᶠᶠ.Cᴮᵀ/Aᶜ

        # Coefficients of the cubic polynomial
        za2 = Pᶜᵒᵉᶠᶠ.Cᴮᵀ*(1 - Rᴮᴬ) + Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*(1-Rᶜᴬ)
        za1 = Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*Pᶜᵒᵉᶠᶠ.Cᴮᵀ*(1 - Rᴮᴬ - Rᶜᴬ) + Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂*(1 - (Rᶜᴬ+Rᶜᴬ))
        za0 = Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂*Pᶜᵒᵉᶠᶠ.Cᴮᵀ*(1 - Rᴮᴬ - (Rᶜᴬ+Rᶜᴬ))

        # Taylor expansion around the minimum 
        #discriminant of the quadratic equation 
        #for the minimum close to the root
        zd = za2*za2 - 3*za1 

        if zd > 0
            if za2 < 0
                zhmin = (-za2 + sqrt(zd))/3
            else
                zhmin = -za1/(za2 + sqrt(zd))
            end

            return zhmin + sqrt(-(za0 + zhmin*(za1 + zhmin*(za2 + zhmin)))/sqrt(zd))
        else
            return 1e-7
        end
    end
end

"""
     FboundsAᵀₙₕ₂ₒ(
        Cᵀ, Pᵀ, Siᵀ, NH₄ᵀ=0, H₂Sᵀ=0, Pᶜᵒᵉᶠᶠ
     )
Calculate the lower and upper bounds of the "non-water-selfionization"
 contributions to total alkalinity.
"""
@inline function FboundsAᵀₙₕ₂ₒ( 
    Cᵀ, Pᵀ, Siᵀ, NH₄ᵀ, H₂Sᵀ, Pᶜᵒᵉᶠᶠ
)
# greatest lower bound (infimum)
    Aᵀₗₒ = -Pᵀ - Pᶜᵒᵉᶠᶠ.Cˢᴼ⁴  - Pᶜᵒᵉᶠᶠ.Cᶠᵀ

# least upper bound (supremum)
    Aᵀₕᵢ = Cᵀ + Cᵀ + Pᶜᵒᵉᶠᶠ.Cᴮᵀ +
                Pᵀ + Pᵀ + Siᵀ +
              NH₄ᵀ + H₂Sᵀ

    return Aᵀₗₒ, Aᵀₕᵢ
end

"""
    function FAᵀ(Cᵀ, Aᵀ, Bᵀ, Pᵀ, Siᵀ,  SO₄ᵀ, Fᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
    
Evaluate the rational function form of the total alkalinity-pH equation
"""
@inline function FAᵀ(
    Cᵀ, Aᵀ, Pᵀ, Siᵀ, NH₄ᵀ, H₂Sᵀ, H⁺, Pᶜᵒᵉᶠᶠ
    )
    
    return FACᵀ(Cᵀ, H⁺, Pᶜᵒᵉᶠᶠ) + 
           FAPᵀ(Pᵀ, H⁺, Pᶜᵒᵉᶠᶠ) + 
           FASiᵀ(Siᵀ, H⁺, Pᶜᵒᵉᶠᶠ) +
           FANH₄ᵀ(NH₄ᵀ, H⁺, Pᶜᵒᵉᶠᶠ) + 
           FAH₂Sᵀ(H₂Sᵀ, H⁺, Pᶜᵒᵉᶠᶠ) +
           FABᵀ(H⁺, Pᶜᵒᵉᶠᶠ) +  
           FASO₄ᵀ(H⁺, Pᶜᵒᵉᶠᶠ) + 
           FAFᵀ(H⁺, Pᶜᵒᵉᶠᶠ) +
           FAH₂O(H⁺, Pᶜᵒᵉᶠᶠ) - 
           Aᵀ,
           F∂A∂Cᵀ(Cᵀ, H⁺, Pᶜᵒᵉᶠᶠ) + 
           F∂A∂Pᵀ(Pᵀ, H⁺, Pᶜᵒᵉᶠᶠ) + 
           F∂A∂Siᵀ(Siᵀ, H⁺, Pᶜᵒᵉᶠᶠ) +
           F∂A∂NH₄ᵀ(NH₄ᵀ, H⁺, Pᶜᵒᵉᶠᶠ) + 
           F∂A∂H₂Sᵀ(H₂Sᵀ, H⁺, Pᶜᵒᵉᶠᶠ) + 
           F∂A∂Bᵀ(H⁺, Pᶜᵒᵉᶠᶠ) + 
           F∂A∂SO₄ᵀ(H⁺, Pᶜᵒᵉᶠᶠ) + 
           F∂A∂Fᵀ(H⁺, Pᶜᵒᵉᶠᶠ) +
           F∂A∂H₂O(H⁺, Pᶜᵒᵉᶠᶠ)
end

"""
    function FACᵀ(Cᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function FACᵀ(Cᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
    # H2CO3 - HCO3 - CO3 : n=2, m=0
    return Cᵀ * (( 2 * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ * 
                       Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ + 
                  H⁺ * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁
                 )/(
                       Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ *
                       Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ + 
                  H⁺*( Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ + H⁺ )
                 ))
end
    
"""
    F∂A∂Cᵀ(Cᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function F∂A∂Cᵀ(Cᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
    # H2CO3 - HCO3 - CO3 : n=2, m=0
    return - Cᵀ * (
                    ( Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ * 
                      Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ *
                      Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ + 
                      H⁺ * ( 4 * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ *
                                 Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ +
                            H⁺ * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁)
                    )/(
                       Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ *
                       Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ + 
                       H⁺*( Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ + H⁺ )
                      )^2
                  ) 
end
    
"""
    function FABᵀ(H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function FABᵀ(H⁺, Pᶜᵒᵉᶠᶠ)
    # B(OH)3 - B(OH)4 : n=1, m=0
    return Pᶜᵒᵉᶠᶠ.Cᴮᵀ * (Pᶜᵒᵉᶠᶠ.Cᵇₖ₁/(Pᶜᵒᵉᶠᶠ.Cᵇₖ₁ + H⁺))   
end

"""
    function F∂A∂Bᵀ(H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function F∂A∂Bᵀ(H⁺, Pᶜᵒᵉᶠᶠ)
    # B(OH)3 - B(OH)4 : n=1, m=0
    return - Pᶜᵒᵉᶠᶠ.Cᴮᵀ * ( 
                            Pᶜᵒᵉᶠᶠ.Cᵇₖ₁
                          )/(
                            Pᶜᵒᵉᶠᶠ.Cᵇₖ₁ + H⁺
                          )^2  
end

"""
    function FAPᵀ(Pᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function FAPᵀ(Pᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
    # H3PO4 - H2PO4 - HPO4 - PO4 : n=3, m=1
    return Pᵀ * 3 * (Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ * 
                     Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂ * 
                     Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₃ + 
               H⁺ * ( 2 * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ * 
                          Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂ + 
                     H⁺ * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁
                    )
                )/(
                 Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ * 
                 Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂ * 
                 Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₃ + 
                 H⁺ * (Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ * 
                       Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂ + 
                       H⁺ * (Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ + H⁺)
                      )
                - 1 ) 
end

"""
    function F∂A∂Pᵀ(Pᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function F∂A∂Pᵀ(Pᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
    # H3PO4 - H2PO4 - HPO4 - PO4 : n=3, m=1
    return - Pᵀ * ((Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ * 
                    Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂ * 
                    Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ * 
                    Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂ *
                    Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₃ + 
                    H⁺ * ( 4 * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ *
                               Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ * 
                               Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂ * 
                               Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₃ + 
                      H⁺ * ( 9 * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ * 
                                 Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂ * 
                                 Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₃ +
                                 Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ *
                                 Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ * 
                                 Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂ +
                        H⁺ * ( 4 * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ * 
                                   Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂ +
                          H⁺ * ( Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ ))))         
                   )/(
                    Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ * 
                    Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂ * 
                    Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₃ + 
                    H⁺ * (Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ * 
                          Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂ + 
                          H⁺ * (Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ + H⁺)
                         )
                   )^2
                  ) 
end

"""
    function FASiᵀ(Siᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function FASiᵀ(Siᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
    # H4SiO4 - H3SiO4 : n=1, m=0
    return Siᵀ * (Pᶜᵒᵉᶠᶠ.Cˢⁱᵗₖ₁/(Pᶜᵒᵉᶠᶠ.Cˢⁱᵗₖ₁ + H⁺))
end

"""
    function F∂A∂Siᵀ(Siᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function F∂A∂Siᵀ(Siᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
    # H4SiO4 - H3SiO4 : n=1, m=0
    return - Siᵀ * (Pᶜᵒᵉᶠᶠ.Cˢⁱᵗₖ₁/(Pᶜᵒᵉᶠᶠ.Cˢⁱᵗₖ₁ + H⁺)^2)
end

"""
    function FASO₄ᵀ(H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function FASO₄ᵀ(H⁺, Pᶜᵒᵉᶠᶠ)
    # HSO4 - SO4 : n=1, m=1
    return Pᶜᵒᵉᶠᶠ.Cˢᴼ⁴ * (Pᶜᵒᵉᶠᶠ.Cᴴˢᴼ⁴ₖ₁/(Pᶜᵒᵉᶠᶠ.Cᴴˢᴼ⁴ₖ₁ + H⁺) - 1)
end

"""
    function F∂A∂SO₄ᵀ(H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function F∂A∂SO₄ᵀ(H⁺, Pᶜᵒᵉᶠᶠ)
    # HSO4 - SO4 : n=1, m=1
    return - Pᶜᵒᵉᶠᶠ.Cˢᴼ⁴ * (Pᶜᵒᵉᶠᶠ.Cᴴˢᴼ⁴ₖ₁/(Pᶜᵒᵉᶠᶠ.Cᴴˢᴼ⁴ₖ₁ + H⁺)^2)
end

"""
    function FAFᵀ(H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function FAFᵀ(H⁺, Pᶜᵒᵉᶠᶠ)
    # HF - F : n=1, m=1
    return Pᶜᵒᵉᶠᶠ.Cᶠᵀ * (Pᶜᵒᵉᶠᶠ.Cᴴᶠₖ₁/(Pᶜᵒᵉᶠᶠ.Cᴴᶠₖ₁ + H⁺) - 1)
end

"""
    function F∂A∂Fᵀ(H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function F∂A∂Fᵀ(H⁺, Pᶜᵒᵉᶠᶠ)
    # HF - F : n=1, m=1
    return - Pᶜᵒᵉᶠᶠ.Cᶠᵀ * (Pᶜᵒᵉᶠᶠ.Cᴴᶠₖ₁/(Pᶜᵒᵉᶠᶠ.Cᴴᶠₖ₁ + H⁺)^2)
end

"""
    function FANH₄ᵀ(NH₄ᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function FANH₄ᵀ(NH₄ᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
    # NH4 - NH3 : n=1, m=0
    return NH₄ᵀ * (Pᶜᵒᵉᶠᶠ.Cᴺᴴ⁴ₖ₁/(Pᶜᵒᵉᶠᶠ.Cᴺᴴ⁴ₖ₁ + H⁺))
end

"""
    function F∂A∂NH₄ᵀ(NH₄ᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function F∂A∂NH₄ᵀ(NH₄ᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
    # NH4 - NH3 : n=1, m=0
    return - NH₄ᵀ * (Pᶜᵒᵉᶠᶠ.Cᴺᴴ⁴ₖ₁/(Pᶜᵒᵉᶠᶠ.Cᴺᴴ⁴ₖ₁ + H⁺)^2)
end

"""
    function FAH₂Sᵀ(H₂Sᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function FAH₂Sᵀ(H₂Sᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
    # H2S - HS : n=1, m=0
    return H₂Sᵀ * (Pᶜᵒᵉᶠᶠ.Cᴴ²ˢₖ₁/(Pᶜᵒᵉᶠᶠ.Cᴴ²ˢₖ₁ + H⁺))
end

"""
    function F∂A∂H₂Sᵀ(H₂Sᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function F∂A∂H₂Sᵀ(H₂Sᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
    # H2S - HS : n=1, m=0
    return - H₂Sᵀ * (Pᶜᵒᵉᶠᶠ.Cᴴ²ˢₖ₁/(Pᶜᵒᵉᶠᶠ.Cᴴ²ˢₖ₁ + H⁺)^2)
end

"""
    function FAH₂O(H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function FAH₂O(H⁺, Pᶜᵒᵉᶠᶠ)
    # H2O - OH
    return Pᶜᵒᵉᶠᶠ.Cᴴ²ᴼₖ₁/H⁺ -H⁺/Pᶜᵒᵉᶠᶠ.H⁺ₜoverH⁺₃
end

"""
    function F∂A∂H₂O(H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function F∂A∂H₂O(H⁺, Pᶜᵒᵉᶠᶠ)
    # H2O - OH
    return - Pᶜᵒᵉᶠᶠ.Cᴴ²ᴼₖ₁/H⁺^2 - 1/Pᶜᵒᵉᶠᶠ.H⁺ₜoverH⁺₃
end

end # module
# ----------------------------------------------------------------------------------
module AlkalinityCorrectionCarbonSolver
export AlkalinityCorrectionCarbonSystem,
        CarbonSystem

using ..CarbonSystemSolvers: CarbonSystem, CarbonChemistryCoefficients, FCᵀCO₂ˢᵒˡ, FCᵀCO₃²⁻, FCᵀHCO₃⁻
#include("carbon_chemistry_coefficients.jl")

"""
    AlkalinityCorrectionCarbonSystem(
            Θ       :: FT = 25.0,
            Sᴬ      :: FT = 35.0,
            Δpᵦₐᵣ   :: FT = 0.0,
            Cᵀ      :: FT = 2050.0e-6,
            Aᵀ      :: FT = 2350.0e-6,
            Pᵀ      :: FT = 1.0e-6,
            Siᵀ     :: FT = 15.0e-6,
            pH      :: FT = 8.0,
            pCO₂ᵃᵗᵐ :: FT = 280.0e-6,
            )

Uses the Follows et al (2006) method to solve the distribution of carbon species
"""
@inline function AlkalinityCorrectionCarbonSystem(
        Θᶜ      :: FT = 25.0,
        Sᴬ      :: FT = 35.0,
        Δpᵦₐᵣ   :: FT = 0.0,
        Cᵀ      :: FT = 2050.0e-6,
        Aᵀ      :: FT = 2350.0e-6,
        Pᵀ      :: FT = 1.0e-6,
        Siᵀ     :: FT = 15.0e-6,
        pH      :: FT = 8.0,
        pCO₂ᵃᵗᵐ :: FT = 280.0e-6) where {FT}

    # CarbonChemistryCoefficients are pretty much all in mol/kg, hence the 1e-6 factors for Cᵀ and Aᵀ
    Cᶜᵒᵉᶠᶠ = CarbonChemistryCoefficients(Θᶜ, Sᴬ, Δpᵦₐᵣ)
    
    # Some logic here about choosing coefficient options, particularly Cᵈⁱᶜ 
    Pᶜᵒᵉᶠᶠ = (Cᵈⁱᶜₖ₀ = Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀,
              Cᵈⁱᶜₖ₁ = Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ₗ₀₀,
              Cᵈⁱᶜₖ₂ = Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ₗ₀₀,
              Cᵇₖ₁   = Cᶜᵒᵉᶠᶠ.Cᵇₖ₁,
              Cᴾᴼ⁴ₖ₁ = Cᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁,
              Cᴾᴼ⁴ₖ₂ = Cᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂,
              Cᴾᴼ⁴ₖ₃ = Cᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₃,
              Cˢⁱᵗₖ₁ = Cᶜᵒᵉᶠᶠ.Cˢⁱᵗₖ₁,
              Cᴴˢᴼ⁴ₖ₁= Cᶜᵒᵉᶠᶠ.Cᴴˢᴼ⁴ₖ₁,
              Cᴴᶠₖ₁  = Cᶜᵒᵉᶠᶠ.Cᴴᶠₖ₁,
              Cᴴ²ᴼₖ₁ = Cᶜᵒᵉᶠᶠ.Cᴴ²ᴼₖ₁,
              Cᴮᵀ    = Cᶜᵒᵉᶠᶠ.Cᴮᵀ,
              Cᶠᵀ    = Cᶜᵒᵉᶠᶠ.Cᶠᵀ,
              Cˢᴼ⁴   = Cᶜᵒᵉᶠᶠ.Cˢᴼ⁴,
    )

    # Calculate pH from Aᵀ and Cᵀ and then calculate the rest of the carbon system
    pH     = Fᵖᴴ⁽ᴬᵀ⁺ᶜᵀ⁾(Aᵀ, Cᵀ, Pᵀ, Siᵀ, pH, Pᶜᵒᵉᶠᶠ)
    CO₂ˢᵒˡ = FCᵀCO₂ˢᵒˡ(Cᵀ, pH, Pᶜᵒᵉᶠᶠ)
    HCO₃⁻  = FCᵀHCO₃⁻(Cᵀ, pH, Pᶜᵒᵉᶠᶠ)
    CO₃²⁻  = FCᵀCO₃²⁻(Cᵀ, pH, Pᶜᵒᵉᶠᶠ)
    pCO₂ᵒᶜᵉ= CO₂ˢᵒˡ / Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀ # correct for fugacity of CO₂ in seawater?

    return CarbonSystem(pH, CO₂ˢᵒˡ, HCO₃⁻, CO₃²⁻, Cᵀ, Aᵀ, pCO₂ᵒᶜᵉ, pCO₂ᵃᵗᵐ)
end # end function

"""
    BO₄H₄⁻(pH, Pᶜᵒᵉᶠᶠ) 

Calculate borate (B(OH)₄⁻) contribution to Aᶜ using salinity as a proxy
"""
@inline function BO₄H₄⁻(pH, Pᶜᵒᵉᶠᶠ) 
    # Calculate H⁺ from pH
    H⁺ = 10^-pH

    return Pᶜᵒᵉᶠᶠ.Cᴮᵀ * Pᶜᵒᵉᶠᶠ.Cᵇₖ₁/(H⁺ + Pᶜᵒᵉᶠᶠ.Cᵇₖ₁)
end

"""
    H₃PO₄(Pᵀ, pH, Pᶜᵒᵉᶠᶠ) 

Calculate orthophosphoric acid (H₃PO₄) contribution to Aᶜ
"""
@inline function H₃PO₄(Pᵀ, pH, Pᶜᵒᵉᶠᶠ) 
    # Calculate H⁺ from pH
    H⁺ = 10^-pH

    return (Pᵀ * H⁺^3) / (
                 H⁺^3 +
                 H⁺^2 * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ +
                 H⁺^1 * (Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂) +
                 H⁺^0 * (Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂ * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₃)
         )
end

"""
    H₂PO₄⁻(Pᵀ, pH, Pᶜᵒᵉᶠᶠ)

Calculate the dihydrogen phosphate (H₂PO₄⁻) contribution to Aᶜ
"""
@inline function H₂PO₄⁻(Pᵀ, pH, Pᶜᵒᵉᶠᶠ)
    # Calculate H⁺ from pH
    H⁺ = 10^-pH

    return (Pᵀ * H⁺^2 * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁) / (
                 H⁺^3 +
                 H⁺^2 * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ +
                 H⁺^1 * (Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂) +
                 H⁺^0 * (Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂ * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₃)
         )
end

"""
    HPO₄²⁻(Pᵀ, pH, Pᶜᵒᵉᶠᶠ) 

Calculate the monohydrogen phosphate (HPO₄²⁻) contribution to Aᶜ
"""
@inline function HPO₄²⁻(Pᵀ, pH, Pᶜᵒᵉᶠᶠ) 
    # Calculate H⁺ from pH
    H⁺ = 10^-pH

    return (Pᵀ * H⁺^1 * (Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂)) / (
                 H⁺^3 +
                 H⁺^2 * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ +
                 H⁺^1 * (Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂) +
                 H⁺^0 * (Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂ * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₃)
         )
end

"""
    PO₄³⁻(Pᵀ, pH, Pᶜᵒᵉᶠᶠ) 

Calculate the phosphate (PO₄³⁻) contribution to Aᶜ
"""
@inline function PO₄³⁻(Pᵀ, pH, Pᶜᵒᵉᶠᶠ) 
    # Calculate H⁺ from pH
    H⁺ = 10^-pH

    return (Pᵀ * H⁺^0 * (Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂ * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₃)) / (
                 H⁺^3 +
                 H⁺^2 * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ +
                 H⁺^1 * (Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂) +
                 H⁺^0 * (Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁ * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂ * Pᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₃)
         )
end

"""
    SiO₄H₃⁻(Siᵀ, pH, Pᶜᵒᵉᶠᶠ) 

Calculate the silicate (SiO(OH)₃⁻) contribution to Aᶜ
"""
@inline function SiO₄H₃⁻(Siᵀ, pH, Pᶜᵒᵉᶠᶠ) 
    # Calculate H⁺ from pH
    H⁺ = 10^-pH

    return Siᵀ * Pᶜᵒᵉᶠᶠ.Cˢⁱᵗₖ₁ / (H⁺ + Pᶜᵒᵉᶠᶠ.Cˢⁱᵗₖ₁)
end

"""
    OH⁻(pH, Pᶜᵒᵉᶠᶠ)

Calculate the hydroxide (OH⁻) contribution to Aᶜ
"""
@inline function OH⁻(pH, Pᶜᵒᵉᶠᶠ) 
    # Calculate H⁺ from pH
    H⁺ = 10^-pH

    return Pᶜᵒᵉᶠᶠ.Cᴴ²ᴼₖ₁ / H⁺
end

"""
    H⁺ᶠʳᵉᵉ(pH, Pᶜᵒᵉᶠᶠ)

Calculate the "Free" H⁺ contribution to Aᶜ
"""
@inline function H⁺ᶠʳᵉᵉ(pH, Pᶜᵒᵉᶠᶠ)
    # Calculate H⁺ from pH
    H⁺ = 10^-pH

    return H⁺ * Pᶜᵒᵉᶠᶠ.Cᴴˢᴼ⁴ₖ₁ / (Pᶜᵒᵉᶠᶠ.Cˢᴼ⁴ + Pᶜᵒᵉᶠᶠ.Cᴴˢᴼ⁴ₖ₁)
end

"""
    HSO₄⁻(pH, Pᶜᵒᵉᶠᶠ)

Calculate the hydrogen sulphate (HSO₄⁻) contribution to Aᶜ
"""
@inline function HSO₄⁻(pH, Pᶜᵒᵉᶠᶠ)
    # Calculate H⁺ from pH
    # H⁺ = 10^-pH

    return Pᶜᵒᵉᶠᶠ.Cˢᴼ⁴ * H⁺ᶠʳᵉᵉ(pH, Pᶜᵒᵉᶠᶠ) / (H⁺ᶠʳᵉᵉ(pH, Pᶜᵒᵉᶠᶠ) + Pᶜᵒᵉᶠᶠ.Cᴴˢᴼ⁴ₖ₁)
end
    
"""
    HF(pH, Pᶜᵒᵉᶠᶠ)

Calculate the hydrogen fluoride (HF) contribution to Aᶜ
"""
@inline function HF(pH, Pᶜᵒᵉᶠᶠ)
    # Calculate H⁺ from pH
    H⁺ = 10^-pH

    return Pᶜᵒᵉᶠᶠ.Cᶠᵀ * H⁺ / (H⁺ + Pᶜᵒᵉᶠᶠ.Cᴴᶠₖ₁)
end

"""
    Fᵖᶜᵒ²⁽ᴬᵀ⁺ᶜᵀ⁾(Aᵀ, Cᵀ, pH, Pᶜᵒᵉᶠᶠ)

Solve for ocean pCO₂ given total Alkalinity and DIC

Estimate H⁺ (hydrogen ion conc) using estimate of Aᶜ, carbonate alkalinity
after (Follows et al., 2006)
"""
@inline function Fᵖᴴ⁽ᴬᵀ⁺ᶜᵀ⁾(Aᵀ, Cᵀ, Pᵀ, Siᵀ, pH, Pᶜᵒᵉᶠᶠ)

# Estimate carbonate alkalinity
    Aᶜ = Aᵀ - 
         BO₄H₄⁻(pH, Pᶜᵒᵉᶠᶠ) - 
         OH⁻(pH, Pᶜᵒᵉᶠᶠ) - 
         SiO₄H₃⁻(Siᵀ, pH, Pᶜᵒᵉᶠᶠ) -
         HPO₄²⁻(Pᵀ, pH, Pᶜᵒᵉᶠᶠ) - 
         2 * PO₄³⁻(Pᵀ, pH, Pᶜᵒᵉᶠᶠ) + 
         H₃PO₄(Pᵀ, pH, Pᶜᵒᵉᶠᶠ) + 
         H⁺ᶠʳᵉᵉ(pH, Pᶜᵒᵉᶠᶠ) +
         HSO₄⁻(pH, Pᶜᵒᵉᶠᶠ) + 
         HF(pH, Pᶜᵒᵉᶠᶠ)

# Evaluate better guess of hydrogen ion conc
    H⁺  = 0.5 * ( 
        ((Cᵀ/Aᶜ)-1) * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ + 
        sqrt(
            (1-(Cᵀ/Aᶜ))^2 * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁^2 -
            4*Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ * (1-2*(Cᵀ/Aᶜ))
            ) 
        )

    # Update pH (may want to iterate on this pH?)
    pH = -log10(H⁺)
    return pH 
end # end function

"""
    Fᵖᴴ⁽ᴬᵀ⁺ᵖᶜᵒ²⁾(Aᵀ, pCO₂, Pᵀ, Siᵀ, pH, Pᶜᵒᵉᶠᶠ)

Solve for ocean DIC given total Alkalinity and pCO₂

Estimate H⁺ (hydrogen ion conc) using estimate of Aᶜ, carbonate alkalinity
after (Follows et al., 2006)
"""
@inline function Fᵖᴴ⁽ᴬᵀ⁺ᵖᶜᵒ²⁾(Aᵀ, pCO₂, Pᵀ, Siᵀ, pH, Pᶜᵒᵉᶠᶠ)

# Estimate carbonate alkalinity
    Aᶜ = Aᵀ - 
         BO₄H₄⁻(pH, Pᶜᵒᵉᶠᶠ) - 
         OH⁻(pH, Pᶜᵒᵉᶠᶠ) - 
         SiO₄H₃⁻(Siᵀ, pH, Pᶜᵒᵉᶠᶠ) -
         HPO₄²⁻(Pᵀ, pH, Pᶜᵒᵉᶠᶠ) - 
         2 * PO₄³⁻(Pᵀ, pH, Pᶜᵒᵉᶠᶠ) + 
         H₃PO₄(Pᵀ, pH, Pᶜᵒᵉᶠᶠ) + 
         H⁺ᶠʳᵉᵉ(pH, Pᶜᵒᵉᶠᶠ) +
         HSO₄⁻(pH, Pᶜᵒᵉᶠᶠ) + 
         HF(pH, Pᶜᵒᵉᶠᶠ)

# Evaluate better guess of hydrogen ion conc
#  take account of fugacity in pco2atm in seawater?
    H⁺  = 0.5 * (
                 ((Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀ * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ * pCO₂) / Aᶜ) + 
             sqrt(
                 ((Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀ * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ * pCO₂) / Aᶜ)^2 + 
              8 * (Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀ * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ * pCO₂) / Aᶜ 
                 )
                )

    # Update pH (may want to iterate on this pH?)
    pH = -log10(H⁺)
    return pH 
end # end function
end # module

# ----------------------------------------------------------------------------------
module DirectCubicCarbonSolver
export DirectCubicCarbonSystem,
        CarbonSystem

using RootSolvers
using ..CarbonSystemSolvers: CarbonSystem, CarbonChemistryCoefficients, FCᵀCO₂ˢᵒˡ, FCᵀCO₃²⁻, FCᵀHCO₃⁻
#include("carbon_chemistry_coefficients.jl")

"""
    DirectCubicCarbonSystem(
            Θ       :: FT = 25.0,
            Sᴬ      :: FT = 35.0,
            Δpᵦₐᵣ   :: FT = 0.0,
            Cᵀ      :: FT = 2050.0e-6,
            Aᵀ      :: FT = 2350.0e-6,
            pH      :: FT = 8.0,
            pCO₂ᵃᵗᵐ :: FT = 280.0e-6,
            )

DirectCubicCarbonSolver solves a cubic equation in terms of [H⁺]; 
Not for serious use, but as a placeholder and for testing purposes
"""
@inline function DirectCubicCarbonSystem(
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

    # Calculate pH, pCO2, and carbon species from Aᵀ and Cᵀ
    pH  = Fᵖᴴ⁽ᴬᵀ⁺ᶜᵀ⁾(Aᵀ, Cᵀ, pH, Pᶜᵒᵉᶠᶠ)

    CO₂ˢᵒˡ = FCᵀCO₂ˢᵒˡ(Cᵀ, pH, Pᶜᵒᵉᶠᶠ)
    HCO₃⁻  = FCᵀHCO₃⁻(Cᵀ, pH, Pᶜᵒᵉᶠᶠ)
    CO₃²⁻  = FCᵀCO₃²⁻(Cᵀ, pH, Pᶜᵒᵉᶠᶠ)
    pCO₂ᵒᶜᵉ= CO₂ˢᵒˡ / Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀ # correct for fugacity of CO₂ in seawater?

    return CarbonSystem(pH, CO₂ˢᵒˡ, HCO₃⁻, CO₃²⁻, Cᵀ, Aᵀ, pCO₂ᵒᶜᵉ, pCO₂ᵃᵗᵐ)
end # end function

"""
    Fᵖᴴ⁽ᴬᵀ⁺ᵖᶜᵒ²⁾(Aᵀ, pCO₂ᵃᵗᵐ, pH, Pᶜᵒᵉᶠᶠ)

Solve for DIC given total Alkalinity and pCO₂
"""
@inline function Fᵖᴴ⁽ᴬᵀ⁺ᵖᶜᵒ²⁾(Aᵀ, pCO₂, pH, Pᶜᵒᵉᶠᶠ)
    # Find the real roots of the polynomial using RootSolvers.jl 
    sol = find_zero(  x -> (
        x^3*(Aᵀ) +
        x^2*(
             Pᶜᵒᵉᶠᶠ.Cᵇₖ₁*
             Aᵀ-
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀*
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*
             pCO₂-
             Pᶜᵒᵉᶠᶠ.Cᵇₖ₁*
             Pᶜᵒᵉᶠᶠ.Cᴮᵀ) + 
        x^1*(
            -Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀*
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*
             Pᶜᵒᵉᶠᶠ.Cᵇₖ₁*
             pCO₂-
           2*Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀*
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂*
             pCO₂
             ) +
        X^0*(
          -2*Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀*
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂*
             Pᶜᵒᵉᶠᶠ.Cᵇₖ₁*
             pCO₂
             )        
        ),
        NewtonsMethodAD{Float64}(10^-(pH)),
        #SecantMethod{Float64}(10^-(pH+0.1), 10^-(pH-0.1)),
        CompactSolution());
    
    if sol.converged == true
        H⁺ = sol.root
        # Update pH
        pH = -log10(H⁺)
        return  pH
    else
        error("DirectCubicCarbonSolver did not converge")
        return nothing
    end
end # end function

"""
    Fᵖᴴ⁽ᴬᵀ⁺ᶜᵀ⁾(Aᵀ, Cᵀ, pH, Pᶜᵒᵉᶠᶠ)

Solve for ocean pCO₂ given total Alkalinity and DIC
"""
@inline function Fᵖᴴ⁽ᴬᵀ⁺ᶜᵀ⁾(Aᵀ, Cᵀ, pH, Pᶜᵒᵉᶠᶠ)
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
        H⁺ = sol.root

        # Update pH
        pH = -log10(H⁺)

        return pH 
    else
        error("DirectCubicCarbonSolver did not converge")
        return nothing
    end
end # end function

end # module DirectCubicCarbonSolver

end