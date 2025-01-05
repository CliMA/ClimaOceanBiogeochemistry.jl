module CarbonSystemSolvers
export CarbonCoefficientParameters, 
    CarbonSolverParameters, 
    CarbonSystemParameters, 
    CarbonChemistryCoefficients, 
    CarbonSystem

# definitions of the structs that hold the carbon system parameters and constructor functions
include("carbon_system_parameters.jl")

# functions to calculate the carbon chemistry coefficients
include("carbon_chemistry_coefficients.jl")

struct CarbonSystem{FT<:Real}
    pH      :: FT
    CO₂ˢᵒˡ  :: FT
    HCO₃⁻   :: FT
    CO₃²⁻   :: FT
    Cᵀ      :: FT
    Aᵀ      :: FT
    pCO₂ᵒᶜᵉ :: FT
    pCO₂ᵃᵗᵐ :: FT
    Pᵈⁱᶜₖₛₒₗₐ :: FT
    Pᵈⁱᶜₖₛₒₗₒ :: FT
    Pᵈⁱᶜₖ₀   :: FT
end
"""
    FCᵀCO₂ˢᵒˡ(Cᵀ, pH, Pᶜᵒᵉᶠᶠ)

Calculate the dissolved and hydrated CO₂ concentration in seawater
given the total carbon concentration Cᵀ, pH, and the carbon chemistry coefficients.
"""
@inline function FCᵀCO₂ˢᵒˡ(Cᵀ::Real, pH::Real, Pᶜᵒᵉᶠᶠ) :: Real
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
@inline function FCᵀHCO₃⁻(Cᵀ::Real, pH::Real, Pᶜᵒᵉᶠᶠ) :: Real
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
@inline function FCᵀCO₃²⁻(Cᵀ::Real, pH::Real, Pᶜᵒᵉᶠᶠ) :: Real
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
@inline function FpCO₂CO₂ˢᵒˡ(pCO₂::Real, Pᶜᵒᵉᶠᶠ) :: Real
    # Perhaps take account of fugacity here?
    return Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀ * pCO₂
end

"""
    FpCO₂HCO₃⁻(pCO₂, pH, Pᶜᵒᵉᶠᶠ)

Calculate the bicarbonate ion concentration in seawater
given the pCO₂, pH, and the carbon chemistry coefficients.
"""
@inline function FpCO₂HCO₃⁻(pCO₂::Real, pH::Real, Pᶜᵒᵉᶠᶠ) :: Real
    # Calculate H⁺ from pH
    H⁺ = 10^-pH

    return (Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀ * pCO₂)/H⁺
end

"""
    FpCO₂CO₃²⁻(pCO₂, pH, Pᶜᵒᵉᶠᶠ)

Calculate the carbonate concentration in seawater
given the pCO₂, pH, and the carbon chemistry coefficients.
"""
@inline function FpCO₂CO₃²⁻(pCO₂::Real, pH::Real, Pᶜᵒᵉᶠᶠ) :: Real
    # Calculate H⁺ from pH
    H⁺ = 10^-pH

    return (Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ * Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀ * pCO₂)/(H⁺^2)
end

# ----------------------------------------------------------------------------------
module UniversalRobustCarbonSolver
export UniversalRobustCarbonSystem,
        CarbonSystem

using ..CarbonSystemSolvers: CarbonCoefficientParameters, 
                             CarbonSolverParameters, 
                             CarbonSystemParameters,
                             CarbonChemistryCoefficients, 
                             CarbonSystem, 
                             FCᵀCO₂ˢᵒˡ, 
                             FCᵀCO₃²⁻, 
                             FCᵀHCO₃⁻

"""
    UniversalRobustCarbonSystem(
            pH      :: Real = 8.0,
            pCO₂ᵃᵗᵐ :: Real = 280.0e-6,
            Θ       :: Real = 25.0,
            Sᴬ      :: Real = 35.0,
            Δpᵦₐᵣ   :: Real = 0.0,
            Cᵀ      :: Real = 2050.0e-6,
            Aᵀ      :: Real = 2350.0e-6,
            Pᵀ      :: Real = 1.0e-6,
            Siᵀ     :: Real = 15.0e-6,
            params,
            )

Uses the Munhoven (2013) SolveSAPHE package to solve the distribution of carbon species
"""
@inline function UniversalRobustCarbonSystem(;
        pH      :: Real = 8.0,
        pCO₂ᵃᵗᵐ :: Real = 280.0e-6,
        Θᶜ      :: Real = 25.0,
        Sᴬ      :: Real = 35.0,
        Δpᵦₐᵣ   :: Real = 0.0,
        Cᵀ      :: Real = 2050.0e-6,
        Aᵀ      :: Real = 2350.0e-6,
        Pᵀ      :: Real = 1.0e-6,
        Siᵀ     :: Real = 15.0e-6,
        NH₄ᵀ    :: Real = 0.0,
        H₂Sᵀ    :: Real = 0.0,
        params  :: CarbonSystemParameters = CarbonSystemParameters(),
        )

    # CarbonChemistryCoefficients are pretty much all in mol/kg, hence the 1e-6 factors for Cᵀ and Aᵀ
    Cᶜᵒᵉᶠᶠ = CarbonChemistryCoefficients(params, Θᶜ, Sᴬ, Δpᵦₐᵣ)
    
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
    iter, pH     = Fᵖᴴᵤₙᵢᵣₒ(Aᵀ, Cᵀ, Pᵀ, Siᵀ, NH₄ᵀ, H₂Sᵀ, pH, Pᶜᵒᵉᶠᶠ, params.Sᵒᵖᵗˢ) 
    if iter == params.Sᵒᵖᵗˢ.Iᴴ⁺ₘₐₓ
        error("UniversalRobustCarbonSystem failed to converge")
    end
    CO₂ˢᵒˡ = FCᵀCO₂ˢᵒˡ(Cᵀ, pH, Pᶜᵒᵉᶠᶠ)
    HCO₃⁻  = FCᵀHCO₃⁻(Cᵀ, pH, Pᶜᵒᵉᶠᶠ)
    CO₃²⁻  = FCᵀCO₃²⁻(Cᵀ, pH, Pᶜᵒᵉᶠᶠ)
    pCO₂ᵒᶜᵉ= CO₂ˢᵒˡ / Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀ # correct for fugacity of CO₂ in seawater?

    return CarbonSystem(
        pH, 
        CO₂ˢᵒˡ, 
        HCO₃⁻, 
        CO₃²⁻, 
        Cᵀ, 
        Aᵀ, 
        pCO₂ᵒᶜᵉ, 
        pCO₂ᵃᵗᵐ, 
        Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖₛₒₗₐ, 
        Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖₛₒₗₒ, 
        Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀,
        )
end # end function

adapt_structure( 
    to, cs::CarbonSystem
    ) = CO₂_flux_parameters(
           adapt(to, cs.pH),
           adapt(to, cs.CO₂ˢᵒˡ),
           adapt(to, cs.HCO₃⁻),
           adapt(to, cs.CO₃²⁻),
           adapt(to, cs.Cᵀ),
           adapt(to, cs.Aᵀ),
           adapt(to, cs.pCO₂ᵒᶜᵉ),
           adapt(to, cs.pCO₂ᵃᵗᵐ),
           adapt(to, cs.Pᵈⁱᶜₖₛₒₗₐ),
           adapt(to, cs.Pᵈⁱᶜₖₛₒₗₒ),
	       adapt(to, cs.Pᵈⁱᶜₖ₀),
)

"""
    Fᵖᴴᵤₙᵢᵣₒ(Aᵀ, Cᵀ, Pᵀ, Siᵀ, pH, Pᶜᵒᵉᶠᶠ)

Calculate the pH of seawater given the total alkalinity Aᵀ, total carbon Cᵀ,
total phosphate Pᵀ, total silicate Siᵀ, and the carbon chemistry coefficients.
Uses the SolveSAPHE package (Munhoven et al., 2013), a universal, robust, pH 
solver that converges from any given initial value.
"""
@inline function Fᵖᴴᵤₙᵢᵣₒ(Aᵀ::Real, 
                          Cᵀ::Real, 
                          Pᵀ::Real, 
                          Siᵀ::Real,
                          NH₄ᵀ::Real, 
                          H₂Sᵀ::Real, 
                          pH::Real, 
                          Pᶜᵒᵉᶠᶠ,  
                          Sᵒᵖᵗˢ,
                          ) :: Tuple{Int, Real}
   
    # Initialize some variables
    #Iᴴ⁺                = 0
    Aᵀᵃᵇˢₘᵢₙ           = floatmax(typeof(Aᵀ))
    H⁺ᶠᵃᶜᵗᵒʳ           = 1.0

    # Get a better initial H+ guess or Calculate H⁺ from the given pH
    H⁺ᵢₙᵢ = ifelse(
                pH == 8, 
                FH⁺ᵢₙᵢ(Aᵀ, Cᵀ, Pᶜᵒᵉᶠᶠ), 
                10^-pH,
            )

    # Calculate initial bounds of H+ concentration
    Aᵀₗₒ, Aᵀₕᵢ = FboundsAᵀₙₕ₂ₒ( 
                    Cᵀ, Pᵀ, Siᵀ, NH₄ᵀ, H₂Sᵀ, Pᶜᵒᵉᶠᶠ
                )

    Δₗₒ = (Aᵀ - Aᵀₗₒ)^2 + 4 * Pᶜᵒᵉᶠᶠ.Cᴴ²ᴼₖ₁/Pᶜᵒᵉᶠᶠ.H⁺ₜoverH⁺₃

    H⁺ₘᵢₙ = ifelse(
                Aᵀ ≥ Aᵀₗₒ, 
                2 * Pᶜᵒᵉᶠᶠ.Cᴴ²ᴼₖ₁/( Aᵀ - Aᵀₗₒ + sqrt(Δₗₒ) ),
                Pᶜᵒᵉᶠᶠ.H⁺ₜoverH⁺₃ *( -(Aᵀ - Aᵀₗₒ) + sqrt(Δₗₒ) )/2,
            )

    Δₕᵢ = (Aᵀ - Aᵀₕᵢ)^2 + 4 * Pᶜᵒᵉᶠᶠ.Cᴴ²ᴼₖ₁/Pᶜᵒᵉᶠᶠ.H⁺ₜoverH⁺₃

    H⁺ₘₐₓ = ifelse(
                Aᵀ ≤ Aᵀₕᵢ, 
                Pᶜᵒᵉᶠᶠ.H⁺ₜoverH⁺₃ *( -(Aᵀ - Aᵀₕᵢ) + sqrt(Δₕᵢ) )/2,
                2 * Pᶜᵒᵉᶠᶠ.Cᴴ²ᴼₖ₁/( Aᵀ - Aᵀₕᵢ + sqrt(Δₕᵢ) ),
            )

    # Initial guess for H⁺
    H⁺ = max(min(H⁺ₘₐₓ, H⁺ᵢₙᵢ), H⁺ₘᵢₙ)             
    #H⁺ = sqrt(H⁺ₘₐₓ * H⁺ₘᵢₙ) # Safer(?) than the above line

    ##while abs(H⁺ᶠᵃᶜᵗᵒʳ) > Δₕ₊
    for Iᴴ⁺ in 1:Sᵒᵖᵗˢ.Iᴴ⁺ₘₐₓ
    # Stop iterations once |\delta{[H]}/[H]| < rdel
    # <=> |(H⁺ - H⁺ₚᵣₑ)/H⁺ₚᵣₑ| = |EXP(-Aᵀᵣₐₜ/(∂Aᵀᵣₐₜ∂H⁺*H⁺ₚᵣₑ)) -1| < rdel
    # |EXP(-Aᵀᵣₐₜ/(∂Aᵀᵣₐₜ∂H⁺*H⁺ₚᵣₑ)) -1| ~ |Aᵀᵣₐₜ/(∂Aᵀᵣₐₜ∂H⁺*H⁺ₚᵣₑ)|
    # Alternatively:
    # |\Delta pH| = |Aᵀᵣₐₜ/(∂Aᵀᵣₐₜ∂H⁺*H⁺ₚᵣₑ*LOG(10))|
    #             ~ 1/LOG(10) * |\Delta [H]|/[H]
    #             < 1/LOG(10) * rdel
    # Hence |Aᵀᵣₐₜ/(∂Aᵀᵣₐₜ∂H⁺*H⁺)| < rdel
    # rdel <-- Δₕ₊
     
    ##    if Iᴴ⁺ ≥ Iᴴ⁺ₘₐₓ
    ##        H⁺ = nothing
    ##        break
    ##    end
    ##     
    ##  # Increase the iteration counter
    ##    Iᴴ⁺ += 1

        # remember for next iteration current H⁺ concentration
        H⁺ₚᵣₑ = H⁺

        Aᵀᵣₐₜ, ∂Aᵀᵣₐₜ∂H⁺ = FAᵀ(
              Cᵀ, Aᵀ, Pᵀ, Siᵀ, NH₄ᵀ, H₂Sᵀ, H⁺, Pᶜᵒᵉᶠᶠ
        )

    ##  if Aᵀᵣₐₜ == 0
    ##      break
    ##  end

        # Adapt bracketing interval
        H⁺ₘᵢₙ = ifelse(
                    Aᵀᵣₐₜ > 0, 
                    H⁺ₚᵣₑ, 
                    H⁺ₘᵢₙ
                )
        H⁺ₘₐₓ = ifelse(
                    Aᵀᵣₐₜ < 0, 
                    H⁺ₚᵣₑ, 
                    H⁺ₘₐₓ
                )

        # Calculate the factor for the next iteration
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

        H⁺ = ifelse(
                abs(H⁺ᶠᵃᶜᵗᵒʳ) > Sᵒᵖᵗˢ.H⁺ᵗʰʳᵉˢʰ,
                H⁺ₚᵣₑ * exp( H⁺ᶠᵃᶜᵗᵒʳ ),
                H⁺ₚᵣₑ + (H⁺ᶠᵃᶜᵗᵒʳ * H⁺ₚᵣₑ),
        )

        # Evaluate if the new H⁺ is within the bounds and adjust
        H⁺, H⁺ᶠᵃᶜᵗᵒʳ = ifelse(
            abs(Aᵀᵣₐₜ) ≥ Aᵀᵃᵇˢₘᵢₙ/2,
            # if the function evaluation at the current point is
            # not decreasing faster than with a bisection step (at least linearly)
            # in absolute value take one bisection step on [ph_min, ph_max]
            # ph_new = (ph_min + ph_max)/2d0
            # In terms of [H]_new:
            # [H]_new = 10**(-ph_new)
            #         = 10**(-(ph_min + ph_max)/2d0)
            #         = SQRT(10**(-(ph_min + phmax)))
            #         = SQRT(H⁺ₘₐₓ * H⁺ₘᵢₙ)
            (sqrt(H⁺ₘₐₓ * H⁺ₘᵢₙ), ( H⁺ - H⁺ₚᵣₑ ) / H⁺ₚᵣₑ),
            ifelse(
                  H⁺ < H⁺ₘᵢₙ,
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
                  (sqrt( H⁺ₚᵣₑ * H⁺ₘᵢₙ ), ( H⁺ - H⁺ₚᵣₑ ) / H⁺ₚᵣₑ ),
                  ifelse(
                      H⁺ > H⁺ₘₐₓ,
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
                      (sqrt( H⁺ₚᵣₑ * H⁺ₘₐₓ ), ( H⁺ - H⁺ₚᵣₑ ) / H⁺ₚᵣₑ),
                      (H⁺, H⁺ᶠᵃᶜᵗᵒʳ),
                  ),
            ),
        )

        #if abs(Aᵀᵣₐₜ) ≥ Aᵀᵃᵇˢₘᵢₙ/2
        ## if the function evaluation at the current point is
        ## not decreasing faster than with a bisection step (at least linearly)
        ## in absolute value take one bisection step on [ph_min, ph_max]
        ## ph_new = (ph_min + ph_max)/2d0
        ## In terms of [H]_new:
        ## [H]_new = 10**(-ph_new)
        ##         = 10**(-(ph_min + ph_max)/2d0)
        ##         = SQRT(10**(-(ph_min + phmax)))
        ##         = SQRT(H⁺ₘₐₓ * H⁺ₘᵢₙ)
        #
        #    H⁺        = sqrt(H⁺ₘₐₓ * H⁺ₘᵢₙ)
        #    H⁺ᶠᵃᶜᵗᵒʳ  = ( H⁺ - H⁺ₚᵣₑ ) / H⁺ₚᵣₑ 
        #else
        ## dAᵀᵣₐₜ/dpH = dAᵀᵣₐₜ/d[H] * d[H]/dpH
        ##           = -∂Aᵀᵣₐₜ∂H⁺ * LOG(10) * [H]
        ## \Delta pH = -Aᵀᵣₐₜ/(∂Aᵀᵣₐₜ∂H⁺*d[H]/dpH) = Aᵀᵣₐₜ/(∂Aᵀᵣₐₜ∂H⁺*[H]*LOG(10))
        ## pH_new = pH_old + \deltapH
        ## [H]_new = 10**(-pH_new)
        ##         = 10**(-pH_old - \Delta pH)
        ##         = [H]_old * 10**(-Aᵀᵣₐₜ/(∂Aᵀᵣₐₜ∂H⁺*[H]_old*LOG(10)))
        ##         = [H]_old * EXP(-LOG(10)*Aᵀᵣₐₜ/(∂Aᵀᵣₐₜ∂H⁺*[H]_old*LOG(10)))
        ##         = [H]_old * EXP(-Aᵀᵣₐₜ/(∂Aᵀᵣₐₜ∂H⁺*[H]_old))
        #
        #    H⁺ᶠᵃᶜᵗᵒʳ = -Aᵀᵣₐₜ / ( ∂Aᵀᵣₐₜ∂H⁺ * H⁺ₚᵣₑ )
        #
        #    H⁺ = ifelse(
        #            abs(H⁺ᶠᵃᶜᵗᵒʳ) > H⁺ᵗʰʳᵉˢʰ,
        #            H⁺ₚᵣₑ * exp( H⁺ᶠᵃᶜᵗᵒʳ ),
        #            H⁺ₚᵣₑ + (H⁺ᶠᵃᶜᵗᵒʳ * H⁺ₚᵣₑ),
        #    )
        #
        #    # if [H]_new < [H]_min
        #    # i.e., if ph_new > ph_max then
        #    # take one bisection step on [ph_prev, ph_max]
        #    # ph_new = (ph_prev + ph_max)/2d0
        #    # In terms of [H]_new:
        #    # [H]_new = 10**(-ph_new)
        #    #         = 10**(-(ph_prev + ph_max)/2d0)
        #    #         = SQRT(10**(-(ph_prev + phmax)))
        #    #         = SQRT([H]_old*10**(-ph_max))
        #    #         = SQRT([H]_old * H⁺ₘᵢₙ)
        #    H⁺, H⁺ᶠᵃᶜᵗᵒʳ = ifelse(
        #        H⁺ < H⁺ₘᵢₙ,
        #        (sqrt( H⁺ₚᵣₑ * H⁺ₘᵢₙ ), ( H⁺ - H⁺ₚᵣₑ ) / H⁺ₚᵣₑ ),
        #        (H⁺, H⁺ᶠᵃᶜᵗᵒʳ)
        #    )
        #
        #    # if [H]_new > [H]_max
        #    # i.e., if ph_new < ph_min, then
        #    # take one bisection step on [ph_min, ph_prev]
        #    # ph_new = (ph_prev + ph_min)/2d0
        #    # In terms of [H]_new:
        #    # [H]_new = 10**(-ph_new)
        #    #         = 10**(-(ph_prev + ph_min)/2d0)
        #    #         = SQRT(10**(-(ph_prev + ph_min)))
        #    #         = SQRT([H]_old*10**(-ph_min))
        #    #         = SQRT([H]_old * zhmax)
        #    H⁺, H⁺ᶠᵃᶜᵗᵒʳ = ifelse(
        #        H⁺ > H⁺ₘₐₓ,
        #        (sqrt( H⁺ₚᵣₑ * H⁺ₘₐₓ ), ( H⁺ - H⁺ₚᵣₑ ) / H⁺ₚᵣₑ),
        #        (H⁺, H⁺ᶠᵃᶜᵗᵒʳ),
        #    )
        #    end
        #end

        if abs(H⁺ᶠᵃᶜᵗᵒʳ) < Sᵒᵖᵗˢ.Δₕ₊
            # H⁺ has converged to the desired accuracy so begin exiting
            # This is similar to what is done in RootSolvers.jl
            Aᵀᵃᵇˢₘᵢₙ = min( abs(Aᵀᵣₐₜ), Aᵀᵃᵇˢₘᵢₙ)
            
            Aᵀᵣₐₜ, ∂Aᵀᵣₐₜ∂H⁺ = ifelse(
                H⁺ > 0,
                FAᵀ(
                    Cᵀ, Aᵀ, Pᵀ, Siᵀ, NH₄ᵀ, H₂Sᵀ, H⁺, Pᶜᵒᵉᶠᶠ
                ),
                (nothing, nothing)
            )
            return convert(AbstractFloat,Iᴴ⁺), -log10(H⁺)
        end

        Aᵀᵃᵇˢₘᵢₙ = min( abs(Aᵀᵣₐₜ), Aᵀᵃᵇˢₘᵢₙ)
    end # end while loop

    Aᵀᵣₐₜ, ∂Aᵀᵣₐₜ∂H⁺ = ifelse(
        H⁺ > 0,
        FAᵀ(
            Cᵀ, Aᵀ, Pᵀ, Siᵀ, NH₄ᵀ, H₂Sᵀ, H⁺, Pᶜᵒᵉᶠᶠ
        ),
        (nothing, nothing)
    )
    return Sᵒᵖᵗˢ.Iᴴ⁺ₘₐₓ, -log10(H⁺)
end

"""
    FH⁺ᵢₙᵢ(Aᵀ, Cᵀ, H⁺, Pᶜᵒᵉᶠᶠ)

Calculates the root for the 2nd order approximation of the
 Cᵀ-Bᵀ-Aᶜ equation for H⁺ (reformulated as a cubic polynomial)
 around the local minimum, if it exists.

 Returns * 1e-03 if Aᶜ <= 0
         * 1e-10 if Aᶜ >= 2*Cᵀ + Bᵀ
         * 1e-07 if 0 < Aᶜ < 2*Cᵀ + Bᵀ
            and the 2nd order approximation does not have a solution

"""
@inline function FH⁺ᵢₙᵢ(Aᶜ::Real, Cᵀ::Real, Pᶜᵒᵉᶠᶠ) :: Real
    # Calculate the coefficients of the cubic polynomial
    Rᶜᴬ = Cᵀ/Aᶜ
    Rᴮᴬ = Pᶜᵒᵉᶠᶠ.Cᴮᵀ/Aᶜ

    # Coefficients of the cubic polynomial
    a₀ = Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂*Pᶜᵒᵉᶠᶠ.Cᴮᵀ*(1 - Rᴮᴬ - (Rᶜᴬ+Rᶜᴬ))
    a₁ = Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*Pᶜᵒᵉᶠᶠ.Cᴮᵀ*(1 - Rᴮᴬ - Rᶜᴬ) + Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂*(1 - (Rᶜᴬ+Rᶜᴬ))
    a₂ = Pᶜᵒᵉᶠᶠ.Cᴮᵀ*(1 - Rᴮᴬ) + Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*(1-Rᶜᴬ)

    # Taylor expansion around the minimum 
    #discriminant of the quadratic equation 
    #for the minimum close to the root
    d = a₂*a₂ - 3*a₁ 

    Hₘᵢₙ = ifelse(
        a₂ < 0,
        (-a₂ + sqrt(d))/3,
        -a₁/(a₂ + sqrt(d))
    )
    
    # Determine a good value for initial H⁺ concentration
    H⁺ᵢₙᵢ = ifelse(
                Aᶜ <= 0,
                1e-3,
                ifelse(
                    Aᶜ >= (2*Cᵀ + Pᶜᵒᵉᶠᶠ.Cᴮᵀ),
                    1e-10,
                    ifelse(
                        d > 0,
                        Hₘᵢₙ + sqrt(-(a₀ + Hₘᵢₙ*(a₁ + Hₘᵢₙ*(a₂ + Hₘᵢₙ)))/sqrt(d)),
                        1e-7,
                    ),
                ),
            )
    return H⁺ᵢₙᵢ
end

"""
     FboundsAᵀₙₕ₂ₒ(
        Cᵀ, Pᵀ, Siᵀ, NH₄ᵀ=0, H₂Sᵀ=0, Pᶜᵒᵉᶠᶠ
     )
Calculate the lower and upper bounds of the "non-water-selfionization"
 contributions to total alkalinity.
"""
@inline function FboundsAᵀₙₕ₂ₒ( 
    Cᵀ::Real, Pᵀ::Real, Siᵀ::Real, NH₄ᵀ::Real, H₂Sᵀ::Real, Pᶜᵒᵉᶠᶠ
) :: Tuple{Real, Real}
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
    Cᵀ::Real, Aᵀ::Real, Pᵀ::Real, Siᵀ::Real, NH₄ᵀ::Real, H₂Sᵀ::Real, H⁺::Real, Pᶜᵒᵉᶠᶠ
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
@inline function FACᵀ(Cᵀ::Real, H⁺::Real, Pᶜᵒᵉᶠᶠ) :: Real
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
@inline function F∂A∂Cᵀ(Cᵀ::Real, H⁺::Real, Pᶜᵒᵉᶠᶠ) :: Real
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
@inline function FABᵀ(H⁺::Real, Pᶜᵒᵉᶠᶠ) :: Real
    # B(OH)3 - B(OH)4 : n=1, m=0
    return Pᶜᵒᵉᶠᶠ.Cᴮᵀ * (Pᶜᵒᵉᶠᶠ.Cᵇₖ₁/(Pᶜᵒᵉᶠᶠ.Cᵇₖ₁ + H⁺))   
end

"""
    function F∂A∂Bᵀ(H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function F∂A∂Bᵀ(H⁺::Real, Pᶜᵒᵉᶠᶠ) :: Real
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
@inline function FAPᵀ(Pᵀ::Real, H⁺::Real, Pᶜᵒᵉᶠᶠ) :: Real
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
@inline function F∂A∂Pᵀ(Pᵀ::Real, H⁺::Real, Pᶜᵒᵉᶠᶠ) :: Real
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
@inline function FASiᵀ(Siᵀ::Real, H⁺::Real, Pᶜᵒᵉᶠᶠ) :: Real
    # H4SiO4 - H3SiO4 : n=1, m=0
    return Siᵀ * (Pᶜᵒᵉᶠᶠ.Cˢⁱᵗₖ₁/(Pᶜᵒᵉᶠᶠ.Cˢⁱᵗₖ₁ + H⁺))
end

"""
    function F∂A∂Siᵀ(Siᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function F∂A∂Siᵀ(Siᵀ::Real, H⁺::Real, Pᶜᵒᵉᶠᶠ) :: Real
    # H4SiO4 - H3SiO4 : n=1, m=0
    return - Siᵀ * (Pᶜᵒᵉᶠᶠ.Cˢⁱᵗₖ₁/(Pᶜᵒᵉᶠᶠ.Cˢⁱᵗₖ₁ + H⁺)^2)
end

"""
    function FASO₄ᵀ(H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function FASO₄ᵀ(H⁺::Real, Pᶜᵒᵉᶠᶠ) :: Real
    # HSO4 - SO4 : n=1, m=1
    return Pᶜᵒᵉᶠᶠ.Cˢᴼ⁴ * (Pᶜᵒᵉᶠᶠ.Cᴴˢᴼ⁴ₖ₁/(Pᶜᵒᵉᶠᶠ.Cᴴˢᴼ⁴ₖ₁ + H⁺) - 1)
end

"""
    function F∂A∂SO₄ᵀ(H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function F∂A∂SO₄ᵀ(H⁺::Real, Pᶜᵒᵉᶠᶠ) :: Real
    # HSO4 - SO4 : n=1, m=1
    return - Pᶜᵒᵉᶠᶠ.Cˢᴼ⁴ * (Pᶜᵒᵉᶠᶠ.Cᴴˢᴼ⁴ₖ₁/(Pᶜᵒᵉᶠᶠ.Cᴴˢᴼ⁴ₖ₁ + H⁺)^2)
end

"""
    function FAFᵀ(H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function FAFᵀ(H⁺::Real, Pᶜᵒᵉᶠᶠ) :: Real
    # HF - F : n=1, m=1
    return Pᶜᵒᵉᶠᶠ.Cᶠᵀ * (Pᶜᵒᵉᶠᶠ.Cᴴᶠₖ₁/(Pᶜᵒᵉᶠᶠ.Cᴴᶠₖ₁ + H⁺) - 1)
end

"""
    function F∂A∂Fᵀ(H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function F∂A∂Fᵀ(H⁺::Real, Pᶜᵒᵉᶠᶠ) :: Real
    # HF - F : n=1, m=1
    return - Pᶜᵒᵉᶠᶠ.Cᶠᵀ * (Pᶜᵒᵉᶠᶠ.Cᴴᶠₖ₁/(Pᶜᵒᵉᶠᶠ.Cᴴᶠₖ₁ + H⁺)^2)
end

"""
    function FANH₄ᵀ(NH₄ᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function FANH₄ᵀ(NH₄ᵀ::Real, H⁺::Real, Pᶜᵒᵉᶠᶠ) :: Real
    # NH4 - NH3 : n=1, m=0
    return NH₄ᵀ * (Pᶜᵒᵉᶠᶠ.Cᴺᴴ⁴ₖ₁/(Pᶜᵒᵉᶠᶠ.Cᴺᴴ⁴ₖ₁ + H⁺))
end

"""
    function F∂A∂NH₄ᵀ(NH₄ᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function F∂A∂NH₄ᵀ(NH₄ᵀ::Real, H⁺::Real, Pᶜᵒᵉᶠᶠ) :: Real
    # NH4 - NH3 : n=1, m=0
    return - NH₄ᵀ * (Pᶜᵒᵉᶠᶠ.Cᴺᴴ⁴ₖ₁/(Pᶜᵒᵉᶠᶠ.Cᴺᴴ⁴ₖ₁ + H⁺)^2)
end

"""
    function FAH₂Sᵀ(H₂Sᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function FAH₂Sᵀ(H₂Sᵀ::Real, H⁺::Real, Pᶜᵒᵉᶠᶠ) :: Real
    # H2S - HS : n=1, m=0
    return H₂Sᵀ * (Pᶜᵒᵉᶠᶠ.Cᴴ²ˢₖ₁/(Pᶜᵒᵉᶠᶠ.Cᴴ²ˢₖ₁ + H⁺))
end

"""
    function F∂A∂H₂Sᵀ(H₂Sᵀ, H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function F∂A∂H₂Sᵀ(H₂Sᵀ::Real, H⁺::Real, Pᶜᵒᵉᶠᶠ) :: Real
    # H2S - HS : n=1, m=0
    return - H₂Sᵀ * (Pᶜᵒᵉᶠᶠ.Cᴴ²ˢₖ₁/(Pᶜᵒᵉᶠᶠ.Cᴴ²ˢₖ₁ + H⁺)^2)
end

"""
    function FAH₂O(H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function FAH₂O(H⁺::Real, Pᶜᵒᵉᶠᶠ) :: Real
    # H2O - OH
    return Pᶜᵒᵉᶠᶠ.Cᴴ²ᴼₖ₁/H⁺ -H⁺/Pᶜᵒᵉᶠᶠ.H⁺ₜoverH⁺₃
end

"""
    function F∂A∂H₂O(H⁺, Pᶜᵒᵉᶠᶠ)
"""
@inline function F∂A∂H₂O(H⁺::Real, Pᶜᵒᵉᶠᶠ) :: Real
    # H2O - OH
    return - Pᶜᵒᵉᶠᶠ.Cᴴ²ᴼₖ₁/H⁺^2 - 1/Pᶜᵒᵉᶠᶠ.H⁺ₜoverH⁺₃
end

end # module
# ----------------------------------------------------------------------------------
module AlkalinityCorrectionCarbonSolver
export AlkalinityCorrectionCarbonSystem,
        CarbonSystem

using ..CarbonSystemSolvers: CarbonCoefficientParameters, 
                             CarbonSolverParameters, 
                             CarbonSystemParameters,
                             CarbonChemistryCoefficients, 
                             CarbonSystem, 
                             FCᵀCO₂ˢᵒˡ, 
                             FCᵀCO₃²⁻, 
                             FCᵀHCO₃⁻

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
@inline function AlkalinityCorrectionCarbonSystem(;
        Θᶜ      :: Real = 25.0,
        Sᴬ      :: Real = 35.0,
        Δpᵦₐᵣ   :: Real = 0.0,
        Cᵀ      :: Real = 2050.0e-6,
        Aᵀ      :: Real = 2350.0e-6,
        Pᵀ      :: Real = 1.0e-6,
        Siᵀ     :: Real = 15.0e-6,
        pH      :: Real = 8.0,
        pCO₂ᵃᵗᵐ :: Real = 280.0e-6,
        params :: CarbonSystemParameters = CarbonSystemParameters()
        )
    
    # CarbonChemistryCoefficients are pretty much all in mol/kg, hence the 1e-6 factors for Cᵀ and Aᵀ
    Cᶜᵒᵉᶠᶠ = CarbonChemistryCoefficients(params, Θᶜ, Sᴬ, Δpᵦₐᵣ)
    
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

    return CarbonSystem(
        pH, 
        CO₂ˢᵒˡ, 
        HCO₃⁻, 
        CO₃²⁻, 
        Cᵀ, 
        Aᵀ, 
        pCO₂ᵒᶜᵉ, 
        pCO₂ᵃᵗᵐ, 
        Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖₛₒₗₐ, 
        Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖₛₒₗₒ, 
        Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀,
        )
end # end function

adapt_structure( 
    to, cs::CarbonSystem
    ) = CO₂_flux_parameters(
           adapt(to, cs.pH),
           adapt(to, cs.CO₂ˢᵒˡ),
           adapt(to, cs.HCO₃⁻),
           adapt(to, cs.CO₃²⁻),
           adapt(to, cs.Cᵀ),
           adapt(to, cs.Aᵀ),
           adapt(to, cs.pCO₂ᵒᶜᵉ),
           adapt(to, cs.pCO₂ᵃᵗᵐ),
           adapt(to, cs.Pᵈⁱᶜₖₛₒₗₐ),
           adapt(to, cs.Pᵈⁱᶜₖₛₒₗₒ),
	       adapt(to, cs.Pᵈⁱᶜₖ₀),
)

"""
    BO₄H₄⁻(pH, Pᶜᵒᵉᶠᶠ) 

Calculate borate (B(OH)₄⁻) contribution to Aᶜ using salinity as a proxy
"""
@inline function BO₄H₄⁻(pH::Real, Pᶜᵒᵉᶠᶠ) :: Real
    # Calculate H⁺ from pH
    H⁺ = 10^-pH

    return Pᶜᵒᵉᶠᶠ.Cᴮᵀ * Pᶜᵒᵉᶠᶠ.Cᵇₖ₁/(H⁺ + Pᶜᵒᵉᶠᶠ.Cᵇₖ₁)
end

"""
    H₃PO₄(Pᵀ, pH, Pᶜᵒᵉᶠᶠ) 

Calculate orthophosphoric acid (H₃PO₄) contribution to Aᶜ
"""
@inline function H₃PO₄(Pᵀ::Real, pH::Real, Pᶜᵒᵉᶠᶠ) :: Real
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
@inline function H₂PO₄⁻(Pᵀ::Real, pH::Real, Pᶜᵒᵉᶠᶠ) :: Real
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
@inline function HPO₄²⁻(Pᵀ::Real, pH::Real, Pᶜᵒᵉᶠᶠ) :: Real
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
@inline function PO₄³⁻(Pᵀ::Real, pH::Real, Pᶜᵒᵉᶠᶠ) :: Real
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
@inline function SiO₄H₃⁻(Siᵀ::Real, pH::Real, Pᶜᵒᵉᶠᶠ) :: Real
    # Calculate H⁺ from pH
    H⁺ = 10^-pH

    return Siᵀ * Pᶜᵒᵉᶠᶠ.Cˢⁱᵗₖ₁ / (H⁺ + Pᶜᵒᵉᶠᶠ.Cˢⁱᵗₖ₁)
end

"""
    OH⁻(pH, Pᶜᵒᵉᶠᶠ)

Calculate the hydroxide (OH⁻) contribution to Aᶜ
"""
@inline function OH⁻(pH::Real, Pᶜᵒᵉᶠᶠ) :: Real
    # Calculate H⁺ from pH
    H⁺ = 10^-pH

    return Pᶜᵒᵉᶠᶠ.Cᴴ²ᴼₖ₁ / H⁺
end

"""
    H⁺ᶠʳᵉᵉ(pH, Pᶜᵒᵉᶠᶠ)

Calculate the "Free" H⁺ contribution to Aᶜ
"""
@inline function H⁺ᶠʳᵉᵉ(pH::Real, Pᶜᵒᵉᶠᶠ) :: Real
    # Calculate H⁺ from pH
    H⁺ = 10^-pH

    return H⁺ * Pᶜᵒᵉᶠᶠ.Cᴴˢᴼ⁴ₖ₁ / (Pᶜᵒᵉᶠᶠ.Cˢᴼ⁴ + Pᶜᵒᵉᶠᶠ.Cᴴˢᴼ⁴ₖ₁)
end

"""
    HSO₄⁻(pH, Pᶜᵒᵉᶠᶠ)

Calculate the hydrogen sulphate (HSO₄⁻) contribution to Aᶜ
"""
@inline function HSO₄⁻(pH::Real, Pᶜᵒᵉᶠᶠ) :: Real
    # Calculate H⁺ from pH
    # H⁺ = 10^-pH

    return Pᶜᵒᵉᶠᶠ.Cˢᴼ⁴ * H⁺ᶠʳᵉᵉ(pH, Pᶜᵒᵉᶠᶠ) / (H⁺ᶠʳᵉᵉ(pH, Pᶜᵒᵉᶠᶠ) + Pᶜᵒᵉᶠᶠ.Cᴴˢᴼ⁴ₖ₁)
end
    
"""
    HF(pH, Pᶜᵒᵉᶠᶠ)

Calculate the hydrogen fluoride (HF) contribution to Aᶜ
"""
@inline function HF(pH::Real, Pᶜᵒᵉᶠᶠ) :: Real
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
@inline function Fᵖᴴ⁽ᴬᵀ⁺ᶜᵀ⁾(Aᵀ::Real, Cᵀ::Real, Pᵀ::Real, Siᵀ::Real, pH::Real, Pᶜᵒᵉᶠᶠ) :: Real

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
@inline function Fᵖᴴ⁽ᴬᵀ⁺ᵖᶜᵒ²⁾(Aᵀ::Real, pCO₂::Real, Pᵀ::Real, Siᵀ::Real, pH::Real, Pᶜᵒᵉᶠᶠ) :: Real

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

using ..CarbonSystemSolvers: CarbonCoefficientParameters, 
                             CarbonSolverParameters, 
                             CarbonSystemParameters,
                             CarbonChemistryCoefficients, 
                             CarbonSystem, 
                             FCᵀCO₂ˢᵒˡ, 
                             FCᵀCO₃²⁻, 
                             FCᵀHCO₃⁻
using RootSolvers

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
@inline function DirectCubicCarbonSystem(;
        Θᶜ      :: Real = 25.0,
        Sᴬ      :: Real = 35.0,
        Δpᵦₐᵣ   :: Real = 0.0,
        Cᵀ      :: Real = 2050.0e-6,
        Aᵀ      :: Real = 2350.0e-6,
        pH      :: Real = 8.0,
        pCO₂ᵃᵗᵐ :: Real = 280.0e-6,
        params :: CarbonSystemParameters = CarbonSystemParameters(),
        )

    # CarbonChemistryCoefficients are pretty much all in mol/kg, hence the 1e-6 factors for Cᵀ and Aᵀ
    Cᶜᵒᵉᶠᶠ = CarbonChemistryCoefficients(params, Θᶜ, Sᴬ, Δpᵦₐᵣ)
    
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

    return CarbonSystem(
        pH, 
        CO₂ˢᵒˡ, 
        HCO₃⁻, 
        CO₃²⁻, 
        Cᵀ, 
        Aᵀ, 
        pCO₂ᵒᶜᵉ, 
        pCO₂ᵃᵗᵐ, 
        Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖₛₒₗₐ, 
        Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖₛₒₗₒ, 
        Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀,
        )
end # end function

adapt_structure( 
    to, cs::CarbonSystem
    ) = CO₂_flux_parameters(
           adapt(to, cs.pH),
           adapt(to, cs.CO₂ˢᵒˡ),
           adapt(to, cs.HCO₃⁻),
           adapt(to, cs.CO₃²⁻),
           adapt(to, cs.Cᵀ),
           adapt(to, cs.Aᵀ),
           adapt(to, cs.pCO₂ᵒᶜᵉ),
           adapt(to, cs.pCO₂ᵃᵗᵐ),
           adapt(to, cs.Pᵈⁱᶜₖₛₒₗₐ),
           adapt(to, cs.Pᵈⁱᶜₖₛₒₗₒ),
	       adapt(to, cs.Pᵈⁱᶜₖ₀),
)

"""
    Fᵖᴴ⁽ᴬᵀ⁺ᵖᶜᵒ²⁾(Aᵀ, pCO₂ᵃᵗᵐ, pH, Pᶜᵒᵉᶠᶠ)

Solve for DIC given total Alkalinity and pCO₂
"""
@inline function Fᵖᴴ⁽ᴬᵀ⁺ᵖᶜᵒ²⁾(Aᵀ::Real, pCO₂::Real, pH::Real, Pᶜᵒᵉᶠᶠ) :: Real
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
        x^0*(
          -2*Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀*
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁*
             Pᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂*
             Pᶜᵒᵉᶠᶠ.Cᵇₖ₁*
             pCO₂
             ),        
        ),
        NewtonsMethodAD{Float64}(10^-(pH)),
        #SecantMethod{Float64}(10^-(pH+0.1), 10^-(pH-0.1)),
        CompactSolution());
    
    
    H⁺ = ifelse(
        sol.converged == true, 
        sol.root, 
        0,
        #error("DirectCubicCarbonSolver did not converge"),
    ) 
    return -log10(H⁺)
end # end function

"""
    Fᵖᴴ⁽ᴬᵀ⁺ᶜᵀ⁾(Aᵀ, Cᵀ, pH, Pᶜᵒᵉᶠᶠ)

Solve for ocean pCO₂ given total Alkalinity and DIC
"""
@inline function Fᵖᴴ⁽ᴬᵀ⁺ᶜᵀ⁾(Aᵀ::Real, Cᵀ::Real, pH::Real, Pᶜᵒᵉᶠᶠ) :: Real
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

    H⁺ = ifelse(
        sol.converged, 
        sol.root, 
        0,
        #error("DirectCubicCarbonSolver did not converge"),
    ) 
    return -log10(H⁺)
end # end function

end # module DirectCubicCarbonSolver

end
