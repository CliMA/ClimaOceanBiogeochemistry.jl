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
module AlkalinityCorrectionCarbonSolver
export AlkalinityCorrectionCarbonSystem

using ..CarbonSystemSolvers: CarbonSystem, FCᵀCO₂ˢᵒˡ, FCᵀCO₃²⁻, FCᵀHCO₃⁻
include("carbon_chemistry_coefficients.jl")

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
export DirectCubicCarbonSystem

using RootSolvers
using ..CarbonSystemSolvers: CarbonSystem, FCᵀCO₂ˢᵒˡ, FCᵀCO₃²⁻, FCᵀHCO₃⁻
include("carbon_chemistry_coefficients.jl")

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

# ----------------------------------------------------------------------------------

# using .CarbonSystemSolvers.DirectCubicCarbonSolver
# using .CarbonSystemSolvers.AlkalinityCorrectionCarbonSolver
# 
# include("carbon_chemistry_coefficients.jl")
# 
# ## This should go in the testing suite, eventually.
# Θᶜ      = 25.0
# Sᴬ      = 35.0
# Δpᵦₐᵣ   = 0.0
# Cᵀ      = 2050e-6 # umol/kg to mol/kg
# Aᵀ      = 2350e-6 # umol/kg to mol/kg
# pCO₂ᵃᵗᵐ = 280e-6  # uatm to atm
# pH      = 8.0
# FT = Float64
# 
# Cᶜᵒᵉᶠᶠ = CarbonChemistryCoefficients(Θᶜ, Sᴬ, Δpᵦₐᵣ)
# """
# #Check the values of calculated constants
# println("Cᵈⁱᶜₖ₀ = ",       Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀     )
# println("Cᵈⁱᶜₖ₁ᵣ₉₃ = ",    Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ᵣ₉₃  )
# println("Cᵈⁱᶜₖ₂ᵣ₉₃ = ",    Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ᵣ₉₃  )
# println("Cᵈⁱᶜₖ₁ₘ₉₅ = ",    Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ₘ₉₅  )
# println("Cᵈⁱᶜₖ₂ₘ₉₅ = ",    Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ₘ₉₅  )
# println("Cᵈⁱᶜₖ₁ₗ₀₀ = ",     Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ₗ₀₀  )
# println("Cᵈⁱᶜₖ₂ₗ₀₀ = ",     Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ₗ₀₀  )
# println("Cᵇₖ₁ = ",         Cᶜᵒᵉᶠᶠ.Cᵇₖ₁       )
# println("Cᴴ²ᴼₖ₁ = ",       Cᶜᵒᵉᶠᶠ.Cᴴ²ᴼₖ₁     )
# println("Cᴾᴼ⁴ₖ₁ = ",       Cᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁     )
# println("Cᴾᴼ⁴ₖ₂ = ",       Cᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂     )
# println("Cᴾᴼ⁴ₖ₃ = ",       Cᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₃     )
# println("Cˢⁱᵗₖ₁ = ",       Cᶜᵒᵉᶠᶠ.Cˢⁱᵗₖ₁     )
# println("Cᴴ²ˢₖ₁ = ",       Cᶜᵒᵉᶠᶠ.Cᴴ²ˢₖ₁     )
# println("Cᴺᴴ⁴ₖ₁ = ",       Cᶜᵒᵉᶠᶠ.Cᴺᴴ⁴ₖ₁     )
# println("Cᴴᶠᵦ₁ = ",        Cᶜᵒᵉᶠᶠ.Cᴴᶠᵦ₁     )
# println("Cᴴᶠₖ₁ = ",        Cᶜᵒᵉᶠᶠ.Cᴴᶠₖ₁      )
# println("Cᴴˢᴼ⁴ₖ₁ = ",      Cᶜᵒᵉᶠᶠ.Cᴴˢᴼ⁴ₖ₁    )
# println("Cᶜᵃˡᶜⁱᵗᵉₛₚ = ",   Cᶜᵒᵉᶠᶠ.Cᶜᵃˡᶜⁱᵗᵉₛₚ  )
# println("Cᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ = ", Cᶜᵒᵉᶠᶠ.Cᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ)
# println("Cᴮᵀ = " ,        Cᶜᵒᵉᶠᶠ.Cᴮᵀ         )
# println("Cᶠᵀ = " ,        Cᶜᵒᵉᶠᶠ.Cᶠᵀ         )
# println("Cᶜᵃ = ",         Cᶜᵒᵉᶠᶠ.Cᶜᵃ         )
# println("Cˢᴼ⁴ = ",        Cᶜᵒᵉᶠᶠ.Cˢᴼ⁴        )
# """
# 
# @assert round(log(Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₀), digits = 4)     == -3.5617 # Handbook (2007)
# @assert round(log(Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ᵣ₉₃), digits = 4)  == -13.4847 # Handbook (1994)
# @assert round(log(Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ᵣ₉₃), digits = 4)  == -20.5504 # Handbook (1994)
# #@assert Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ₘ₉₅  ==
# #@assert Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ₘ₉₅  ==
# @assert round(log10(Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₁ₗ₀₀), digits = 4)   == -5.8472 # Handbook (2007)
# @assert round(log10(Cᶜᵒᵉᶠᶠ.Cᵈⁱᶜₖ₂ₗ₀₀), digits = 4)   == -8.9660 # Handbook (2007)
# @assert round(log(Cᶜᵒᵉᶠᶠ.Cᵇₖ₁), digits = 4)       == -19.7964 # Handbook (2007)
# @assert round(log(Cᶜᵒᵉᶠᶠ.Cᴴ²ᴼₖ₁*Cᶜᵒᵉᶠᶠ.H⁺ₛoverH⁺ₜ)-0.015, digits = 3)     == -30.434 # Handbook (2007)
# @assert round(log(Cᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₁*Cᶜᵒᵉᶠᶠ.H⁺ₛoverH⁺ₜ)-0.015, digits = 2)     == -3.71 # Handbook (2007)
# @assert round(log(Cᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₂*Cᶜᵒᵉᶠᶠ.H⁺ₛoverH⁺ₜ)-0.015, digits = 3)     == -13.727 # Handbook (2007)
# @assert round(log(Cᶜᵒᵉᶠᶠ.Cᴾᴼ⁴ₖ₃*Cᶜᵒᵉᶠᶠ.H⁺ₛoverH⁺ₜ)-0.015, digits = 2)     == -20.24 # Handbook (2007)
# @assert round(log(Cᶜᵒᵉᶠᶠ.Cˢⁱᵗₖ₁*Cᶜᵒᵉᶠᶠ.H⁺ₛoverH⁺ₜ)-0.015, digits = 2)     == -21.61 # Handbook (2007)
# @assert round(-log10(Cᶜᵒᵉᶠᶠ.Cᴴ²ˢₖ₁*Cᶜᵒᵉᶠᶠ.H⁺ₛoverH⁺ₜ), digits = 2)     == 6.51 # Lewis and Wallace (1998)
# @assert round(-log10(Cᶜᵒᵉᶠᶠ.Cᴺᴴ⁴ₖ₁*Cᶜᵒᵉᶠᶠ.H⁺ₛoverH⁺ₜ), digits = 2)     == 9.26 # Lewis and Wallace (1998)
# #@assert Cᶜᵒᵉᶠᶠ.Cᴴᶠᵦ₁      ==
# @assert round(log(Cᶜᵒᵉᶠᶠ.Cᴴᶠₖ₁), digits = 2)       == -6.09 # Handbook (2007)
# @assert round(log(Cᶜᵒᵉᶠᶠ.Cᴴˢᴼ⁴ₖ₁), digits = 2)     == -2.30 # Handbook (2007)
# #@assert Cᶜᵒᵉᶠᶠ.Cᶜᵃˡᶜⁱᵗᵉₛₚ  ==
# #@assert Cᶜᵒᵉᶠᶠ.Cᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ==
# #@assert Cᶜᵒᵉᶠᶠ.Cᴮᵀ        ==
# #@assert Cᶜᵒᵉᶠᶠ.Cᶠᵀ        ==
# #@assert Cᶜᵒᵉᶠᶠ.Cᶜᵃ        ==
# #@assert Cᶜᵒᵉᶠᶠ.Cˢᴼ⁴       ==
# println("Testing DirectCubicCarbonSolver for pCO2:")
# (; pH, CO₂ˢᵒˡ, HCO₃⁻, CO₃²⁻, Cᵀ, Aᵀ, pCO₂ᵒᶜᵉ, pCO₂ᵃᵗᵐ) = 
# DirectCubicCarbonSystem(
#         Θᶜ, Sᴬ, Δpᵦₐᵣ, Cᵀ, Aᵀ, pH, pCO₂ᵃᵗᵐ,
#         )
# 
# println("Cᵀ = ", Cᵀ * 1e6  )
# println("Aᵀ = ", Aᵀ * 1e6  )
# println("pH = " , pH       )
# println("pCO₂ᵃᵗᵐ = ", pCO₂ᵃᵗᵐ * 1e6)
# println("pCO₂ᵒᶜᵉ = ", pCO₂ᵒᶜᵉ * 1e6)
# println("")
# 
# Pᵀ = 0.5e-6  # umol/kg to mol/kg
# Siᵀ = 7.5e-6 # umol/kg to mol/kg
# 
# println("Testing AlkalinityCorrectionCarbonSolver (Follows et al., 2006) for pCO2:")
# (; pH, CO₂ˢᵒˡ, HCO₃⁻, CO₃²⁻, Cᵀ, Aᵀ, pCO₂ᵒᶜᵉ, pCO₂ᵃᵗᵐ) = 
# AlkalinityCorrectionCarbonSystem(
#         Θᶜ, Sᴬ, Δpᵦₐᵣ, Cᵀ, Aᵀ, Pᵀ, Siᵀ, pH, pCO₂ᵃᵗᵐ,
#         )
# 
# println("Cᵀ = ", Cᵀ * 1e6  )
# println("Aᵀ = ", Aᵀ * 1e6  )
# println("Pᵀ = ", Pᵀ * 1e6  )
# println("Siᵀ = ", Siᵀ * 1e6  )
# println("pH = " , pH       )
# println("pCO₂ᵃᵗᵐ = ", pCO₂ᵃᵗᵐ * 1e6)
# println("pCO₂ᵒᶜᵉ = ", pCO₂ᵒᶜᵉ * 1e6)
# println("NB, you wouldn't expect to get exactly the same results here \nbecause of the extra terms in the calcite alkalinity.")
end