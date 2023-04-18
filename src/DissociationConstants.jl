module DissociationConstants

export DissociationCoefficients

struct DissociationCoefficients{FT}
    """
    DissociationCoefficients{FT}
    Returns dissociation coefficients necessary to solve for the distribution of carbonate species.
    
    "Temperature in degrees Celsius"
    Θ     :: FT
    "Absolute Salinity in g/kg"
    Sᴬ    :: FT
    "Applied Pressure in bars"
    Δpᵦₐᵣ :: FT
    """
    function DissociationCoefficients(
            Θᶜ    :: FT = 25.0,
            Sᴬ    :: FT = 35.0,
            Δpᵦₐᵣ :: FT = 0.0,
            ) where {FT}

# Need a conversion from Absolute Salinity, Sᴬ, to Practical Salinity, Sᴾ
    Sᵖ = Sᴬ
# Sᴾ = Sᴬ / (1.0 - 0.35 * Sᴬ / 35.0) ??
# What about converting temperature from Conservative Temperature to potential temperature?
    Θᴷ = Θᶜ + Θᴷ_0ᵒC
# Also need to convert from Absolute Pressure, Pᴬ, to Applied Pressure in bars, the pressure relative to (1x) atm

    return (Cᵈⁱᶜₖ₀    = Fᵈⁱᶜₖ₀(Θᴷ, Sᵖ, Pᵈⁱᶜₖ₀),
            Cᵈⁱᶜₖ₁ᵣ₉₃ = Fᵈⁱᶜₖ₁ᵣ₉₃(Θᴷ, Sᵖ, Δpᵦₐᵣ, Pᵈⁱᶜₖ₁ᵣ₉₃),
            Cᵈⁱᶜₖ₂ᵣ₉₃ = Fᵈⁱᶜₖ₂ᵣ₉₃(Θᴷ, Sᵖ, Δpᵦₐᵣ, Pᵈⁱᶜₖ₂ᵣ₉₃),
            Cᵈⁱᶜₖ₁ₘ₉₅ = Fᵈⁱᶜₖ₁ₘ₉₅(Θᴷ, Sᵖ, Δpᵦₐᵣ, Pᵈⁱᶜₖ₁ₘ₉₅)/H⁺ₛoverH⁺ₜ(Θᴷ, Sᵖ, Δpᵦₐᵣ),
            Cᵈⁱᶜₖ₂ₘ₉₅ = Fᵈⁱᶜₖ₂ₘ₉₅(Θᴷ, Sᵖ, Δpᵦₐᵣ, Pᵈⁱᶜₖ₂ₘ₉₅)/H⁺ₛoverH⁺ₜ(Θᴷ, Sᵖ, Δpᵦₐᵣ),
            Cᵈⁱᶜₖ₁ₗ₀₀  = Fᵈⁱᶜₖ₁ₗ₀₀(Θᴷ, Sᵖ, Δpᵦₐᵣ, Pᵈⁱᶜₖ₁ₗ₀₀),
            Cᵈⁱᶜₖ₂ₗ₀₀  = Fᵈⁱᶜₖ₂ₗ₀₀(Θᴷ, Sᵖ, Δpᵦₐᵣ, Pᵈⁱᶜₖ₂ₗ₀₀),
            Cᵇₖ₁      = Fᵇₖ₁(Θᴷ, Sᵖ, Δpᵦₐᵣ, Pᴮₖ₁),
            Cᴴ²ᴼₖ₁    = Fᴴ²ᴼₖ₁(Θᴷ, Sᵖ, Δpᵦₐᵣ, Pᴴ²ᴼₖ₁)/H⁺ₛoverH⁺ₜ(Θᴷ, Sᵖ, Δpᵦₐᵣ),
            Cᴾᴼ⁴ₖ₁    = Fᴾᴼ⁴ₖ₁(Θᴷ, Sᵖ, Δpᵦₐᵣ, Pᴾᴼ⁴ₖ₁)/H⁺ₛoverH⁺ₜ(Θᴷ, Sᵖ, Δpᵦₐᵣ),
            Cᴾᴼ⁴ₖ₂    = Fᴾᴼ⁴ₖ₂(Θᴷ, Sᵖ, Δpᵦₐᵣ, Pᴾᴼ⁴ₖ₂)/H⁺ₛoverH⁺ₜ(Θᴷ, Sᵖ, Δpᵦₐᵣ),
            Cᴾᴼ⁴ₖ₃    = Fᴾᴼ⁴ₖ₃(Θᴷ, Sᵖ, Δpᵦₐᵣ, Pᴾᴼ⁴ₖ₃)/H⁺ₛoverH⁺ₜ(Θᴷ, Sᵖ, Δpᵦₐᵣ),
            Cˢⁱᵗₖ₁    = Fˢⁱᵗₖ₁(Θᴷ, Sᵖ, Pˢⁱᵗₖ₁)/H⁺ₛoverH⁺ₜ(Θᴷ, Sᵖ, Δpᵦₐᵣ),
            Cᴴ²ˢₖ₁    = Fᴴ²ˢₖ₁(Θᴷ, Sᵖ, Δpᵦₐᵣ, Pᴴ²ˢₖ₁)/H⁺ₛoverH⁺ₜ(Θᴷ, Sᵖ, Δpᵦₐᵣ),
            Cᴺᴴ⁴ₖ₁    = Fᴺᴴ⁴ₖ₁(Θᴷ, Sᵖ, Δpᵦₐᵣ, Pᴺᴴ⁴ₖ₁)/H⁺ₛoverH⁺ₜ(Θᴷ, Sᵖ, Δpᵦₐᵣ),
            Cᴴᶠᵦ₁     = Fᴴᶠᵦ₁(Θᴷ, Sᵖ, Δpᵦₐᵣ, Pᴴᶠᵦ₁),
            Cᴴᶠₖ₁     = Fᴴᶠₖ₁(Θᴷ, Sᵖ, Δpᵦₐᵣ, Pᴴᶠₖ₁),
            Cᴴˢᴼ⁴ₖ₁   = Fᴴˢᴼ⁴ₖ₁(Θᴷ, Sᵖ, Δpᵦₐᵣ, Pᴴˢᴼ⁴ₖ₁),
            Cᶜᵃˡᶜⁱᵗᵉₛₚ = Fᶜᵃˡᶜⁱᵗᵉₛₚ(Θᴷ, Sᵖ, Δpᵦₐᵣ, Pᶜᵃˡᶜⁱᵗᵉₛₚ),
            Cᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ = Fᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ(Θᴷ, Sᵖ, Δpᵦₐᵣ, Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ),
            Cᴮᵀ         = Bᵀᴼᵀ(Sᵖ),
            Cᶠᵀ         = Fᵀᴼᵀ(Sᵖ),
            Cᶜᵃ         = Caᵀᴼᵀ(Sᵖ),
            Cˢᴼ⁴        = SO₄ᵀᴼᵀ(Sᵖ),
            CH⁺ₛoverH⁺ₜ  = H⁺ₛoverH⁺ₜ(Θᴷ, Sᵖ, Δpᵦₐᵣ),
            CH⁺ₜoverH⁺₃ = H⁺ₜoverH⁺₃(Θᴷ, Sᵖ, Δpᵦₐᵣ),
            CH⁺ₛoverH⁺₃ = H⁺ₛoverH⁺₃(Θᴷ, Sᵖ, Δpᵦₐᵣ),
    )
    end # end function
end # end struct

# Gas constant
const gasconst_bar_cm3_o_mol_k = 83.14472 # Handbook (2007)

# 0 degrees centigrade in Kelvin
const Θᴷ_0ᵒC = 273.15 # Handbook (2007)

@inline Θᴷ₁₀₀(Θᴷ)    = Θᴷ/100
@inline Sᴾ⁰⁵(Sᵖ)     = sqrt(Sᵖ)
@inline ΘᵒC(Θᴷ)      = Θᴷ - Θᴷ_0ᵒC
Base.@kwdef struct Pᵈˢᵖ{FT}
    a₀ :: FT    = - 34.8
end
@inline ΔSᵖ(Sᵖ,Pᵈˢᵖ) = Sᵖ + Pᵈˢᵖ.a₀
@inline Rₜ(Θᴷ)        = gasconst_bar_cm3_o_mol_k * Θᴷ

Base.@kwdef struct Pᴴ²⁰ˢʷ{FT}
    a₀ :: FT    =   1.0
    a₁ :: FT    = - 0.001005
end

@inline 
"""
H₂Oˢʷ(Sᵖ)

Return the mass of pure water in one kg of seawater
of practical salinity, `Sᵖ`.
References: "libthdyct" -- derived by Munhoven (1997) from data by Millero (1982)
            "Handbook (2007)" -- Handbook (2007)
pH scale:   N/A
"""
H₂Oˢʷ(Sᵖ,Pᴴ²⁰ˢʷ) =  Pᴴ²⁰ˢʷ.a₀ + Pᴴ²⁰ˢʷ.a₁ * Sᵖ 
#H₂Oˢʷ(Sᵖ) = 1. - 0.0010049*Sᵖ # libthdyct

Base.@kwdef struct Pᵘˢ{FT}
    a₀ :: FT    =   0.019924
end

@inline 
"""
μₛ(Sᵖ)

Return ionic strength in mol/kg-SW, for given practical salinity, `Sᵖ`.
References: "libthdyct" -- derived by Munhoven (1997) from data by Millero (1982)
            "Handbook (2007)" -- Handbook (2007)
pH scale:   N/A
"""
μₛ(Sᵖ,Pᵘˢ)    = (Pᵘˢ.a₀ * Sᵖ) / H₂Oˢʷ(Sᵖ) # Handbook (2007)
# μₛ(Sᵖ)    = (0.019920 * Sᵖ) / H₂Oˢʷ(Sᵖ)# libthdyct

Base.@kwdef struct Pᴮᵀᴼᵀ{FT}
    a₀ :: FT    =   0.000416
    a₁ :: FT    =   35.0
    a₂ :: FT    =   1.0
end

@inline 
"""
Bᵀᴼᵀ(Sᵖ)

Return total borate concentration in mol/kg-SW given practical salinity, `Sᵖ`.
References: Uppström (1974), cited by  Dickson et al. (2007, chapter 5, p 10)
            Millero (1982) cited in Millero (1995)
"""
Bᵀᴼᵀ(Sᵖ,Pᴮᵀᴼᵀ)  = Pᴮᵀᴼᵀ.a₀ * (Sᵖ / Pᴮᵀᴼᵀ.a₁) / Pᴮᵀᴼᵀ.a₂
#Bᵀᴼᵀ(Sᵖ)  = 0.000232 * (Sᵖ/1.80655)/10.811

Base.@kwdef struct Pᶜᵃᵀᴼᵀ{FT}
    a₀ :: FT    =   0.02127
    a₁ :: FT    =   40.078
    a₂ :: FT    =   1.80655
end

@inline
"""
Caᵀᴼᵀ(Sᵖ)

Return calcium concentration in mol/kg-SW given practical salinity, `Sᵖ`.
References: Culkin and Cox (1966), 
            Culkin (1967), 
            Riley and Tongudai (1967)    
"""
Caᵀᴼᵀ(Sᵖ,Pᶜᵃᵀᴼᵀ) = (Pᶜᵃᵀᴼᵀ.a₀/Pᶜᵃᵀᴼᵀa₁) * (Sᵖ/Pᶜᵃᵀᴼᵀ.a₂)
# Caᵀᴼᵀ(Sᵖ) = 0.010282*(Sᵖ/35.)

Base.@kwdef struct Pᶠᵀᴼᵀ{FT}
    a₀ :: FT    =   0.000068
    a₁ :: FT    =   35.0
end

@inline
"""
Fᵀᴼᵀ(Sᵖ)

Return total fluoride concentration in mol/kg-SW given practical salinity, `Sᵖ`.
References: Culkin (1965) (???)
"""
Fᵀᴼᵀ(Sᵖ,Pᶠᵀᴼᵀ)  = Pᶠᵀᴼᵀ.a₀ * (Sᵖ / Pᶠᵀᴼᵀ.a₁)

Base.@kwdef struct Pˢᴼ⁴ᵀᴼᵀ{FT}
    a₀ :: FT    =   0.1400
    a₁ :: FT    =   96.062
    a₂ :: FT    =   1.80655
end

@inline
"""
SO₄ᵀᴼᵀ(Sᵖ)

Return total sulfate concentration in mol/kg-SW given practical salinity, `Sᵖ`.
References: Morris, A.W. and Riley, J.P. (1966) quoted in Handbook (2007)
"""
SO₄ᵀᴼᵀ(Sᵖ, Pˢᴼ⁴ᵀᴼᵀ) = (Pˢᴼ⁴ᵀᴼᵀ.a₀ / Pˢᴼ⁴ᵀᴼᵀ.a₁ ) * ( Sᵖ / Pˢᴼ⁴ᵀᴼᵀ.a₂)
#SO₄ᵀᴼᵀ(Sᵖ) = 0.028234*(Sᵖ/35.)

Base.@kwdef struct Pᴺᴴ⁴ᵀᴼᵀ{FT}
    a₀ :: FT    =   0.0001
    a₁ :: FT    =   17.031
    a₂ :: FT    =   1.80655
end

@inline
"""
NH₄ᵀᴼᵀ(Sᵖ)

Return total ammonium concentration in mol/kg-SW given practical salinity, `Sᵖ`.
References: Yamamoto (1995)
"""
NH₄ᵀᴼᵀ(Sᵖ,Pᴺᴴ⁴ᵀᴼᵀ) = (Pᴺᴴ⁴ᵀᴼᵀ.a₀ / Pᴺᴴ⁴ᵀᴼᵀ.a₁) * (Sᵖ / Pᴺᴴ⁴ᵀᴼᵀ.a₂)

Base.@kwdef struct Pᴾᴼ⁴ᵀᴼᵀ{FT}
    a₀ :: FT    =   0.0001
    a₁ :: FT    =   94.971
    a₂ :: FT    =   1.80655
end

@inline
"""
PO₄ᵀᴼᵀ(Sᵖ)

Return total phosphate concentration in mol/kg-SW given practical salinity, `Sᵖ`.
References: Millero (1995)
"""
PO₄ᵀᴼᵀ(Sᵖ, Pᴾᴼ⁴ᵀᴼᵀ) = (Pᴾᴼ⁴ᵀᴼᵀ.a₀ / Pᴾᴼ⁴ᵀᴼᵀ.a₁) * (Sᵖ / Pᴾᴼ⁴ᵀᴼᵀ.a₂)

Base.@kwdef struct Pˢⁱᴼ³ᵀᴼᵀ{FT}
    a₀ :: FT    =   0.0001
    a₁ :: FT    =   60.084
    a₂ :: FT    =   1.80655
end

@inline
"""
SiO₃ᵀᴼᵀ(Sᵖ)

Return total silicate concentration in mol/kg-SW given practical salinity, `Sᵖ`.
References: Millero (1995)
"""
SiO₃ᵀᴼᵀ(Sᵖ,Pˢⁱᴼ³ᵀᴼᵀ) = (Pˢⁱᴼ³ᵀᴼᵀ.a₀ / Pˢⁱᴼ³ᵀᴼᵀ.a₁) * (Sᵖ / Pˢⁱᴼ³ᵀᴼᵀ.a₂)

Base.@kwdef struct Pᴴ²ˢᵀᴼᵀ{FT}
    a₀ :: FT    =   0.0001
    a₁ :: FT    =   34.082
    a₂ :: FT    =   1.80655
end

@inline
"""
H₂Sᵀᴼᵀ(Sᵖ)

Return total hydrogen sulfide concentration in mol/kg-SW given practical salinity, `Sᵖ`.
References: Dickson (1990)
"""
H₂Sᵀᴼᵀ(Sᵖ,Pᴴ²ˢᵀᴼᵀ) = (Pᴴ²ˢᵀᴼᵀ.a₀ / Pᴴ²ˢᵀᴼᵀ.a₁) * (Sᵖ / Pᴴ²ˢᵀᴼᵀ.a₂)

Base.@kwdef struct Pᵈⁱᶜₖ₀{FT}
    a₀ :: FT = -60.2409
    a₁ :: FT =  93.4517
    a₂ :: FT =  23.3585
    b₀ :: FT =   0.023517
    b₁ :: FT = - 0.023656
    b₂ :: FT =   0.0047036
end

@inline
"""
Cᵈⁱᶜₖ₀(Θᴷ,Sᵖ,Pᵈⁱᶜₖ₀)

Return hydration constant of CO₂ in (mol/kg-SW)/atm given temperature 
in K, `Θᴷ`, practical salinity, `Sᵖ`, and coefficients, `Pᵈⁱᶜₖ₀`.

CO₂ + H₂O <-> H₂CO₃

References: Weiss (1979)
pH scale  : N/A
Note      : currently no pressure correction
"""
function Fᵈⁱᶜₖ₀(Θᴷ,Sᵖ,Pᵈⁱᶜₖ₀)

    exp( 
        Pᵈⁱᶜₖ₀.a₀ + 
        Pᵈⁱᶜₖ₀.a₁/Θᴷ₁₀₀(Θᴷ) +
        Pᵈⁱᶜₖ₀.a₂*log(Θᴷ₁₀₀(Θᴷ)) +    
        ( 
         Pᵈⁱᶜₖ₀.b₀ + 
         Pᵈⁱᶜₖ₀.b₁*(Θᴷ₁₀₀(Θᴷ)) +
         Pᵈⁱᶜₖ₀.b₂*Θᴷ₁₀₀(Θᴷ)*Θᴷ₁₀₀(Θᴷ)
        )*Sᵖ
       )
end

Base.@kwdef struct Pᵈⁱᶜₖ₁ᵣ₉₃{FT}
    a₀ :: FT =     2.83655
    a₁ :: FT = -2307.1266
    a₂ :: FT = -   1.5529413
    b₀ :: FT = -   0.20760841
    b₁ :: FT = -   4.0484
    b₂ :: FT =     0.08468345
    b₃ :: FT = -   0.00654208
    v₀ :: FT = -  25.5
    v₁ :: FT = -   0.151
    v₂ :: FT =     0.1271
    k₀ :: FT = -   3.08e-3
    k₁ :: FT = -   0.578e-3
    k₂ :: FT =     0.0877e-3
end

@inline
"""
Cᵈⁱᶜₖ₁ᵣ₉₃(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᵈⁱᶜₖ₁ᵣ₉₃)

Return the first dissociation constant of carbonic acid in mol/kg-SW, given
temperature in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, 
`Δpᵦₐᵣ`, and coefficients, `Pᵈⁱᶜₖ₁ᵣ₉₃`.

H₂CO₃ <-> HCO₃⁻ + H⁺

References: Roy et al. (1993) -- also Handbook (1994)
            Millero (1979) pressure correction
pH scale  : Total
Note      : converted here from mol/kg-H2O to mol/kg-SW
"""
function Fᵈⁱᶜₖ₁ᵣ₉₃(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᵈⁱᶜₖ₁ᵣ₉₃)

    ln_kc1_p0 = ( Pᵈⁱᶜₖ₁ᵣ₉₃.a₀ + 
                   Pᵈⁱᶜₖ₁ᵣ₉₃.a₁/Θᴷ + 
                   Pᵈⁱᶜₖ₁ᵣ₉₃.a₂*log(Θᴷ)
                + (
                   Pᵈⁱᶜₖ₁ᵣ₉₃.b₀ + 
                   Pᵈⁱᶜₖ₁ᵣ₉₃.b₁/Θᴷ
                   )*Sᴾ⁰⁵(Sᵖ)
                + Pᵈⁱᶜₖ₁ᵣ₉₃.b₂*Sᵖ
                + Pᵈⁱᶜₖ₁ᵣ₉₃.b₃*Sᵖ*Sᴾ⁰⁵(Sᵖ)
                )

    "Pressure correction for applied pressure /= 0"
    ln_kc1_pp = (-(Pᵈⁱᶜₖ₁ᵣ₉₃.v₀ + 
                   Pᵈⁱᶜₖ₁ᵣ₉₃.v₁*ΔSᵖ(Sᵖ) +
                   Pᵈⁱᶜₖ₁ᵣ₉₃.v₂*ΘᵒC(Θᴷ)) 
                 +(
                     Pᵈⁱᶜₖ₁ᵣ₉₃.k₀ + 
                     Pᵈⁱᶜₖ₁ᵣ₉₃.k₁*ΔSᵖ(Sᵖ) +
                     Pᵈⁱᶜₖ₁ᵣ₉₃.k₂*ΘᵒC(Θᴷ)
                 )*(Δpᵦₐᵣ/2.)
                 )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    exp( ln_kc1_p0 + ln_kc1_pp ) * H₂Oˢʷ(Sᵖ)
end

Base.@kwdef struct Pᵈⁱᶜₖ₂ᵣ₉₃{FT}
    a₀ :: FT = -   9.226508
    a₁ :: FT = -3351.6106
    a₂ :: FT = -   0.2005743
    b₀ :: FT = -   0.106901773
    b₁ :: FT = -  23.9722
    b₂ :: FT =     0.1130822
    b₃ :: FT = -   0.00846934
    v₀ :: FT = -  15.82
    v₁ :: FT =     0.321
    v₂ :: FT = -   0.0219
    k₀ :: FT =     1.13e-3
    k₁ :: FT = -   0.314e-3
    k₂ :: FT = -   0.1475e-3
end

@inline 
"""
Cᵈⁱᶜₖ₂ᵣ₉₃(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᵈⁱᶜₖ₂ᵣ₉₃)

Return the second dissociation constant of carbonic acid in mol/kg-SW, given
temperature in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and
coefficients, `Pᵈⁱᶜₖ₂ᵣ₉₃`.

HCO₃⁻ <-> CO₃²⁻ + H⁺

References: Roy et al. (1993) -- also Handbook (1994)
            Millero (1979) pressure correction
pH scale  : Total
Note      : converted here from mol/kg-H2O to mol/kg-SW
"""
function Fᵈⁱᶜₖ₂ᵣ₉₃(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᵈⁱᶜₖ₂ᵣ₉₃)
    ln_kc2_p0 = ( Pᵈⁱᶜₖ₂ᵣ₉₃.a₀ + 
                   Pᵈⁱᶜₖ₂ᵣ₉₃.a₁/Θᴷ + 
                   Pᵈⁱᶜₖ₂ᵣ₉₃.a₂*log(Θᴷ)
                + (
                   Pᵈⁱᶜₖ₂ᵣ₉₃.b₀ + 
                   Pᵈⁱᶜₖ₂ᵣ₉₃.b₁/Θᴷ
                   )*Sᴾ⁰⁵(Sᵖ)
                +  Pᵈⁱᶜₖ₂ᵣ₉₃.b₂*Sᵖ
                + Pᵈⁱᶜₖ₂ᵣ₉₃.b₃*Sᵖ*Sᴾ⁰⁵(Sᵖ)
                )    

    "Pressure correction for applied pressure /= 0"
    ln_kc2_pp = (-(Pᵈⁱᶜₖ₂ᵣ₉₃.v₀ + 
                   Pᵈⁱᶜₖ₂ᵣ₉₃.v₁*ΔSᵖ(Sᵖ) +
                   Pᵈⁱᶜₖ₂ᵣ₉₃.v₂*ΘᵒC(Θᴷ)) 
                 +(
                     Pᵈⁱᶜₖ₂ᵣ₉₃.k₀ + 
                     Pᵈⁱᶜₖ₂ᵣ₉₃.k₁*ΔSᵖ(Sᵖ) +
                     Pᵈⁱᶜₖ₂ᵣ₉₃.k₂*ΘᵒC(Θᴷ)
                 )*(Δpᵦₐᵣ/2.)
                 )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    exp(ln_kc2_p0 + ln_kc2_pp) * H₂Oˢʷ(Sᵖ)
end

Base.@kwdef struct Pᵈⁱᶜₖ₁ₘ₉₅{FT}
    a₀ :: FT =     2.18867
    a₁ :: FT = -2275.0360
    a₂ :: FT = -   1.468591
    b₀ :: FT = -   0.138681
    b₁ :: FT = -   9.33291
    b₂ :: FT =     0.0726483
    b₃ :: FT = -   0.00574938
    v₀ :: FT = -  25.5
    v₁ :: FT = -   0.151
    v₂ :: FT =     0.1271
    k₀ :: FT = -   3.08e-3
    k₁ :: FT = -   0.578e-3
    k₂ :: FT =     0.0877e-3
end

@inline 
"""
Cᵈⁱᶜₖ₁ₘ₉₅(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᵈⁱᶜₖ₁ₘ₉₅)

Return the first dissociation constant of carbonic acid in mol/kg-SW, given
temperature in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, 
`Δpᵦₐᵣ`, and coefficients, `Pᵈⁱᶜₖ₁ₘ₉₅`.
    
H₂CO₃ <-> HCO₃⁻ + H⁺
    
References: Millero (1995, eq 50 -- ln K1(COM))
             Millero (1982) pressure correction
pH scale:   SWS
"""
function Fᵈⁱᶜₖ₁ₘ₉₅(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᵈⁱᶜₖ₁ₘ₉₅)

    ln_kc1_p0 = ( Pᵈⁱᶜₖ₁ₘ₉₅.a₀ + 
                   Pᵈⁱᶜₖ₁ₘ₉₅.a₁/Θᴷ + 
                   Pᵈⁱᶜₖ₁ₘ₉₅.a₂*log(Θᴷ)
                + (
                   Pᵈⁱᶜₖ₁ₘ₉₅.b₀ + 
                   Pᵈⁱᶜₖ₁ₘ₉₅.b₁/Θᴷ
                   )*Sᴾ⁰⁵(Sᵖ)
                +  Pᵈⁱᶜₖ₁ₘ₉₅.b₂*Sᵖ
                + Pᵈⁱᶜₖ₁ₘ₉₅.b₃*Sᵖ*Sᴾ⁰⁵(Sᵖ)
                )


    "Pressure correction for applied pressure /= 0"
    ln_kc1_pp = (-(Pᵈⁱᶜₖ₁ₘ₉₅.v₀ + 
                   Pᵈⁱᶜₖ₁ₘ₉₅.v₁*ΔSᵖ(Sᵖ) +
                   Pᵈⁱᶜₖ₁ₘ₉₅.v₂*ΘᵒC(Θᴷ)) 
                 +(
                     Pᵈⁱᶜₖ₁ₘ₉₅.k₀ + 
                     Pᵈⁱᶜₖ₁ₘ₉₅.k₁*ΔSᵖ(Sᵖ) +
                     Pᵈⁱᶜₖ₁ₘ₉₅.k₂*ΘᵒC(Θᴷ)
                 )*(Δpᵦₐᵣ/2.)
                 )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    exp( ln_kc1_p0 + ln_kc1_pp )
end

Base.@kwdef struct Pᵈⁱᶜₖ₂ₘ₉₅{FT}
    a₀ :: FT = -   0.84226
    a₁ :: FT = -3741.1288
    a₂ :: FT = -   1.437139
    b₀ :: FT = -   0.128417
    b₁ :: FT = -  24.41239
    b₂ :: FT =     0.1195308
    b₃ :: FT = -   0.00912840
    v₀ :: FT = -  15.82
    v₁ :: FT =     0.321
    v₂ :: FT = -   0.0219
    k₀ :: FT =     1.13e-3
    k₁ :: FT = -   0.314e-3
    k₂ :: FT = -   0.1475e-3
end

@inline
"""
Cᵈⁱᶜₖ₂ₘ₉₅(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᵈⁱᶜₖ₂ₘ₉₅)

Return the second dissociation constant of carbonic acid in mol/kg-SW, given
temperature in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and
coefficients, `Pᵈⁱᶜₖ₂ₘ₉₅`.
    
HCO₃⁻ <-> CO₃²⁻ + H⁺

References: Millero (1995, eq 51 -- ln K2(COM))
            Millero (1979) pressure correction
pH scale:   SWS
"""
function Fᵈⁱᶜₖ₂ₘ₉₅(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᵈⁱᶜₖ₂ₘ₉₅)
    ln_kc2_p0 = (
        Pᵈⁱᶜₖ₂ₘ₉₅.a₀ + 
        Pᵈⁱᶜₖ₂ₘ₉₅.a₁/Θᴷ + 
        Pᵈⁱᶜₖ₂ₘ₉₅.a₂*log(Θᴷ)
        + (
            Pᵈⁱᶜₖ₂ₘ₉₅.b₀ + 
            Pᵈⁱᶜₖ₂ₘ₉₅.b₁/Θᴷ
        )*Sᴾ⁰⁵(Sᵖ) +
        Pᵈⁱᶜₖ₂ₘ₉₅.b₂*Sᵖ +
        Pᵈⁱᶜₖ₂ₘ₉₅.b₃*Sᵖ*Sᴾ⁰⁵(Sᵖ)
    )

    "Pressure correction for applied pressure /= 0"
    ln_kc2_pp = (
        -(Pᵈⁱᶜₖ₂ₘ₉₅.v₀ + 
          Pᵈⁱᶜₖ₂ₘ₉₅.v₁*ΔSᵖ(Sᵖ) +
          Pᵈⁱᶜₖ₂ₘ₉₅.v₂*ΘᵒC(Θᴷ)) 
        +(
            Pᵈⁱᶜₖ₂ₘ₉₅.k₀ + 
            Pᵈⁱᶜₖ₂ₘ₉₅.k₁*ΔSᵖ(Sᵖ) +
            Pᵈⁱᶜₖ₂ₘ₉₅.k₂*ΘᵒC(Θᴷ)
        )*(Δpᵦₐᵣ/2.)
    )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    exp( ln_kc2_p0 + ln_kc2_pp )
end

Base.@kwdef struct Pᵈⁱᶜₖ₁ₗ₀₀{FT}
    a₀ :: FT =    61.2172
    a₁ :: FT = -3633.86
    a₂ :: FT = -   9.67770
    b₀ :: FT =     0.011555
    b₁ :: FT = -   0.0001152
    v₀ :: FT = -  25.5
    v₁ :: FT = -   0.151
    v₂ :: FT =     0.1271
    k₀ :: FT = -   3.08e-3
    k₁ :: FT = -   0.578e-3
    k₂ :: FT =     0.0877e-3
end

@inline
"""
Cᵈⁱᶜₖ₁ₗ₀₀(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᵈⁱᶜₖ₁ₗ₀₀)

Return the first dissociation constant of carbonic acid in mol/kg-SW, given
temperature in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, 
`Δpᵦₐᵣ`, and coefficients, `Pᵈⁱᶜₖ₁ₗ₀₀`.
    
H₂CO₃ <-> HCO₃⁻ + H⁺
    
References: Luecker et al. (2000) -- also Handbook (2007)
            Millero (1979) pressure correction
pH scale:   Total
"""
function Fᵈⁱᶜₖ₁ₗ₀₀(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᵈⁱᶜₖ₁ₗ₀₀)
    log10_kc1_p0 = (
        Pᵈⁱᶜₖ₁ₗ₀₀.a₀ + 
        Pᵈⁱᶜₖ₁ₗ₀₀.a₁/Θᴷ + 
        Pᵈⁱᶜₖ₁ₗ₀₀.a₂*log(Θᴷ)
        + (
            Pᵈⁱᶜₖ₁ₗ₀₀.b₀ + 
            Pᵈⁱᶜₖ₁ₗ₀₀.b₁*Sᵖ
        )*Sᵖ
    )

    "Pressure correction for applied pressure /= 0"
    ln_kc1_pp = (
        -(Pᵈⁱᶜₖ₁ₗ₀₀.v₀ + 
          Pᵈⁱᶜₖ₁ₗ₀₀.v₁*ΔSᵖ(Sᵖ) +
          Pᵈⁱᶜₖ₁ₗ₀₀.v₂*ΘᵒC(Θᴷ)) 
        +(
            Pᵈⁱᶜₖ₁ₗ₀₀.k₀ + 
            Pᵈⁱᶜₖ₁ₗ₀₀.k₁*ΔSᵖ(Sᵖ) +
            Pᵈⁱᶜₖ₁ₗ₀₀.k₂*ΘᵒC(Θᴷ)
        )*(Δpᵦₐᵣ/2.)
    )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    10^(log10_kc1_p0) * exp(ln_kc1_pp)
end

Base.@kwdef struct Pᵈⁱᶜₖ₂ₗ₀₀{FT}
    a₀ :: FT = -  25.9290
    a₁ :: FT = - 471.78
    a₂ :: FT =     3.16967
    b₀ :: FT =     0.01781
    b₁ :: FT = -   0.0001122
    v₀ :: FT = -  15.82
    v₁ :: FT =     0.321
    v₂ :: FT = -   0.0219
    k₀ :: FT =     1.13e-3
    k₁ :: FT = -   0.314e-3
    k₂ :: FT = -   0.1475e-3
end

@inline
"""
Cᵈⁱᶜₖ₂ₗ₀₀(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᵈⁱᶜₖ₂ₗ₀₀)

Return the second dissociation constant of carbonic acid in mol/kg-SW, given
temperature in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and
coefficients, `Pᵈⁱᶜₖ₂ₗ₀₀`.
    
HCO₃⁻ <-> CO₃²⁻ + H⁺

References: Luecker et al. (2000) -- also Handbook (2007)
            Millero (1979) pressure correction
pH scale:   Total
"""
function Fᵈⁱᶜₖ₂ₗ₀₀(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᵈⁱᶜₖ₂ₗ₀₀)
    log10_kc2_p0 = (
        Pᵈⁱᶜₖ₂ₗ₀₀.a₀ + 
        Pᵈⁱᶜₖ₂ₗ₀₀.a₁/Θᴷ + 
        Pᵈⁱᶜₖ₂ₗ₀₀.a₂*log(Θᴷ)
        + (
            Pᵈⁱᶜₖ₂ₗ₀₀.b₀ + 
            Pᵈⁱᶜₖ₂ₗ₀₀.b₁*Sᵖ
        )*Sᵖ
    )

    "Pressure correction for applied pressure /= 0"
    ln_kc2_pp = (
        -(Pᵈⁱᶜₖ₂ₗ₀₀.v₀ + 
          Pᵈⁱᶜₖ₂ₗ₀₀.v₁*ΔSᵖ(Sᵖ) +
          Pᵈⁱᶜₖ₂ₗ₀₀.v₂*ΘᵒC(Θᴷ)) 
        +(
            Pᵈⁱᶜₖ₂ₗ₀₀.k₀ + 
            Pᵈⁱᶜₖ₂ₗ₀₀.k₁*ΔSᵖ(Sᵖ) +
            Pᵈⁱᶜₖ₂ₗ₀₀.k₂*ΘᵒC(Θᴷ)
        )*(Δpᵦₐᵣ/2.)
    )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    10^(log10_kc2_p0) * exp(ln_kc2_pp)
end

Base.@kwdef struct Pᴮₖ₁{FT}
    a₀ :: FT = - 8966.90
    a₁ :: FT = - 2890.53
    a₂ :: FT = -   77.942
    a₃ :: FT =      1.728
    a₄ :: FT = -    0.0996
    b₀ :: FT =    148.0248
    b₁ :: FT =    137.1942
    b₂ :: FT =      1.62142
    c₀ :: FT = -   24.4344
    c₁ :: FT = -   25.085
    c₂ :: FT = -    0.2474
    d₀ :: FT =      0.053105
    v₀ :: FT = -   29.48
    v₁ :: FT =      0.295
    v₂ :: FT =      0.1622
    v₃ :: FT = -    0.002608
    k₀ :: FT = -    2.84e-3
    k₁ :: FT =      0.354e-3
end

@inline
"""
Cᵇₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴮₖ₁)

Return boric acid dissociation constant in mol/kg-SW, given temperature in K,
`Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᴮₖ₁`.

References: Dickson (1990, eq. 23) -- also Handbook (2007, eq. 37)
            Millero (1979) pressure correction
pH scale  : total
"""
function Fᵇₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴮₖ₁)
    ln_kb_p0  = (( Pᴮₖ₁.a₀ +
                Sᴾ⁰⁵(Sᵖ)*( Pᴮₖ₁.a₁ +
                Sᴾ⁰⁵(Sᵖ)*( Pᴮₖ₁.a₂ +
                Sᴾ⁰⁵(Sᵖ)*( Pᴮₖ₁.a₃ + 
                Sᴾ⁰⁵(Sᵖ)*  Pᴮₖ₁.a₄
                    )))) / Θᴷ
                + Pᴮₖ₁.b₀ + 
                Sᴾ⁰⁵(Sᵖ)*( 
                    Pᴮₖ₁.b₁ + 
                    Pᴮₖ₁.b₂*Sᴾ⁰⁵(Sᵖ)
                ) +
                ( Pᴮₖ₁.c₀ + 
                Sᴾ⁰⁵(Sᵖ)*(
                    Pᴮₖ₁.c₁ +
                    Pᴮₖ₁.c₂*Sᴾ⁰⁵(Sᵖ))
                    ) * log(Θᴷ)
                + Pᴮₖ₁.d₀*Sᴾ⁰⁵(Sᵖ)*Θᴷ
                )

    "Pressure correction for applied pressure /= 0"
    ln_kb_pp = (-( Pᴮₖ₁.v₀ + 
                   Pᴮₖ₁.v₁*ΔSᵖ(Sᵖ) +
                   Pᴮₖ₁.v₂*ΘᵒC(Θᴷ) + 
                   Pᴮₖ₁.v₃*ΘᵒC(Θᴷ)*ΘᵒC(Θᴷ)) 
                 +(
                   Pᴮₖ₁.k₀ + 
                   Pᴮₖ₁.k₁*ΔSᵖ(Sᵖ)
                 )*(Δpᵦₐᵣ/2.)
                 )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    exp( ln_kb_p0 + ln_kb_pp )
end

Base.@kwdef struct Pᴴ²ᴼₖ₁{FT}
    a₀ :: FT =    148.9802
    a₁ :: FT = -13847.26
    a₂ :: FT = -   23.6521
    b₀ :: FT = -    5.977
    b₁ :: FT =    118.67
    b₂ :: FT =      1.0495
    c₀ :: FT = -    0.01615
    v₀ :: FT =  -  20.02
    v₁ :: FT =      0.1119
    v₂ :: FT =  -   0.1409E-02
    k₀ :: FT = -    5.13e-3
    k₁ :: FT =      0.0794e-3
end

@inline """
Cᴴ²ᴼₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴ²ᴼₖ₁)

Return dissociation constant of water in (mol/kg-SW)^2, given temperature in K,
`Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᴴ²ᴼₖ₁`.

References: Millero (1995) for value at p_bar = 0
            Millero (pers. comm. 1996) for pressure correction
pH scale  : SWS
"""
function Fᴴ²ᴼₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴ²ᴼₖ₁)
    ln_kw_p0 = (Pᴴ²ᴼₖ₁.a₀ +
                Pᴴ²ᴼₖ₁.a₁/Θᴷ +
                Pᴴ²ᴼₖ₁.a₂*log(Θᴷ) +
                ( Pᴴ²ᴼₖ₁.b₀ +
                Pᴴ²ᴼₖ₁.b₁/Θᴷ +
                Pᴴ²ᴼₖ₁.b₂*log(Θᴷ)
                )*Sᴾ⁰⁵(Sᵖ) +
                Pᴴ²ᴼₖ₁.c₀*Sᵖ
                )

    "Pressure correction for applied pressure /= 0"
    ln_kw_pp = (-( Pᴴ²ᴼₖ₁.v₀ + 
                   Pᴴ²ᴼₖ₁.v₁*ΘᵒC(Θᴷ) + 
                   Pᴴ²ᴼₖ₁.v₂*ΘᵒC(Θᴷ)*ΘᵒC(Θᴷ)) 
                 +(
                   Pᴴ²ᴼₖ₁.k₀ + 
                   Pᴴ²ᴼₖ₁.k₁*ΘᵒC(Θᴷ)
                 )*(Δpᵦₐᵣ/2.)
                 )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    exp( ln_kw_p0 + ln_kw_pp )
end

Base.@kwdef struct Pᴾᴼ⁴ₖ₁{FT}
    a₀ :: FT =   115.54
    a₁ :: FT = -4576.752
    a₂ :: FT = -  18.453
    b₀ :: FT =     0.69171
    b₁ :: FT = - 106.736
    b₂ :: FT = -   0.01844
    b₃ :: FT = -   0.65643
    v₀ :: FT = -  14.51
    v₁ :: FT =     0.1211
    v₂ :: FT = -   0.321E-03
    k₀ :: FT = -   2.67e-3
    k₁ :: FT =     0.0427e-3
end

@inline
"""
Cᴾᴼ⁴ₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴾᴼ⁴ₖ₁)

Return the first dissociation constant of phosphoric acid (H3PO4) in seawater, given 
temperature in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᴾᴼ⁴ₖ₁`.

References: Yao and Millero (1995)
            Millero (1995) for pressure correction
pH scale  : SWS
"""
function Fᴾᴼ⁴ₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴾᴼ⁴ₖ₁)
    ln_kp1_p0 = ( Pᴾᴼ⁴ₖ₁.a₀ +
                  Pᴾᴼ⁴ₖ₁.a₁/Θᴷ +
                  Pᴾᴼ⁴ₖ₁.a₂*log(Θᴷ) +
                 ( 
                  Pᴾᴼ⁴ₖ₁.b₀ +
                  Pᴾᴼ⁴ₖ₁.b₁/Θᴷ
                 )*Sᴾ⁰⁵(Sᵖ) +
                 ( 
                  Pᴾᴼ⁴ₖ₁.b₂ +
                  Pᴾᴼ⁴ₖ₁.b₃/Θᴷ
                 )*Sᵖ )

    "Pressure correction for applied pressure /= 0"
    ln_kp1_pp = (-( Pᴾᴼ⁴ₖ₁.v₀ + 
                    Pᴾᴼ⁴ₖ₁.v₁*ΘᵒC(Θᴷ) +
                    Pᴾᴼ⁴ₖ₁.v₂*ΘᵒC(Θᴷ)*ΘᵒC(Θᴷ)) 
                 +(
                   Pᴾᴼ⁴ₖ₁.k₀ + 
                   Pᴾᴼ⁴ₖ₁.k₁*ΘᵒC(Θᴷ)
                 )*(Δpᵦₐᵣ/2.)
                 )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    exp(ln_kp1_p0 + ln_kp1_pp)
end

Base.@kwdef struct Pᴾᴼ⁴ₖ₂{FT}
    a₀ :: FT =    172.1033
    a₁ :: FT = - 8814.715
    a₂ :: FT = -   27.927
    b₀ :: FT =      1.3566
    b₁ :: FT = -  160.340
    b₂ :: FT = -    0.05778
    b₃ :: FT =      0.37335
    v₀ :: FT = -   23.12
    v₁ :: FT =      0.1758
    v₂ :: FT = -    0.002647
    k₀ :: FT = -    5.15e-3
    k₁ :: FT =      0.09e-3
end

@inline
"""
Cᴾᴼ⁴ₖ₂(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴾᴼ⁴ₖ₂)

Return the second dissociation constant of phosphoric acid (H3PO4) in seawater, given
temperature in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᴾᴼ⁴ₖ₂`.

References: Yao and Millero (1995)
            Millero (1995) for pressure correction
pH scale  : SWS
"""
function Fᴾᴼ⁴ₖ₂(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴾᴼ⁴ₖ₂)
    ln_kp2_p0 = ( Pᴾᴼ⁴ₖ₂.a₀ +
                  Pᴾᴼ⁴ₖ₂.a₁/Θᴷ +
                  Pᴾᴼ⁴ₖ₂.a₂*log(Θᴷ) +
                 ( 
                  Pᴾᴼ⁴ₖ₂.b₀ +
                  Pᴾᴼ⁴ₖ₂.b₁/Θᴷ
                 )*Sᴾ⁰⁵(Sᵖ) +
                 ( 
                  Pᴾᴼ⁴ₖ₂.b₂ +
                  Pᴾᴼ⁴ₖ₂.b₃/Θᴷ
                 )*Sᵖ )

    "Pressure correction for applied pressure /= 0"
    ln_kp2_pp = (-( Pᴾᴼ⁴ₖ₂.v₀ + 
                    Pᴾᴼ⁴ₖ₂.v₁*ΘᵒC(Θᴷ) +
                    Pᴾᴼ⁴ₖ₂.v₂*ΘᵒC(Θᴷ)*ΘᵒC(Θᴷ)) 
                 +(
                   Pᴾᴼ⁴ₖ₂.k₀ + 
                   Pᴾᴼ⁴ₖ₂.k₁*ΘᵒC(Θᴷ)
                 )*(Δpᵦₐᵣ/2.)
                 )*(Δpᵦₐᵣ/Rₜ(Θᴷ))
    
    exp( ln_kp2_p0 + ln_kp2_pp )
end

Base.@kwdef struct Pᴾᴼ⁴ₖ₃{FT}
    a₀ :: FT = -   18.126
    a₁ :: FT = - 3070.75
    a₂ :: FT =      2.81197
    a₃ :: FT =     17.27039
    a₄ :: FT = -    0.09984
    a₅ :: FT = -   44.99486
    v₀ :: FT = -   26.57
    v₁ :: FT =      0.2020
    v₂ :: FT = -    3.042e-3
    k₀ :: FT = -    4.08e-3
    k₁ :: FT =      0.0714e-3
end

@inline
"""
Cᴾᴼ⁴ₖ₃(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴾᴼ⁴ₖ₃)

Return the third dissociation constant of phosphoric acid (H3PO4) in seawater, given 
temperature in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᴾᴼ⁴ₖ₃`.

References: Yao and Millero (1995)
            Millero (1995) for pressure correction
pH scale  : SWS
"""
function Fᴾᴼ⁴ₖ₃(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴾᴼ⁴ₖ₃)
    ln_kp3_p0 = ( Pᴾᴼ⁴ₖ₃.a₀ +
                  Pᴾᴼ⁴ₖ₃.a₁/Θᴷ +
                 ( 
                  Pᴾᴼ⁴ₖ₃.a₂ +
                  Pᴾᴼ⁴ₖ₃.a₃/Θᴷ
                 )*Sᴾ⁰⁵(Sᵖ) +
                 ( 
                  Pᴾᴼ⁴ₖ₃.a₄ +
                  Pᴾᴼ⁴ₖ₃.a₅/Θᴷ
                 )*Sᵖ )

    "Pressure correction for applied pressure /= 0"
    ln_kp3_pp = (-( Pᴾᴼ⁴ₖ₃.v₀ + 
                    Pᴾᴼ⁴ₖ₃.v₁*ΘᵒC(Θᴷ) +
                    Pᴾᴼ⁴ₖ₃.v₂*ΘᵒC(Θᴷ)*ΘᵒC(Θᴷ)) 
                 +(
                   Pᴾᴼ⁴ₖ₃.k₀ + 
                   Pᴾᴼ⁴ₖ₃.k₁*ΘᵒC(Θᴷ)
                 )*(Δpᵦₐᵣ/2.)
                 )*(Δpᵦₐᵣ/Rₜ(Θᴷ))
    
    exp( ln_kp3_p0 + ln_kp3_pp )
end

Base.@kwdef struct Pˢⁱᵗₖ₁{FT}
    a₀ :: FT =    117.40
    a₁ :: FT = - 8904.2
    a₂ :: FT = -   19.334
    b₀ :: FT =      3.5913
    b₁ :: FT = -  458.79
    b₂ :: FT = -    1.5998
    b₃ :: FT =    188.74
    c₀ :: FT =      0.07871
    c₁ :: FT = -   12.1652
end

@inline
"""
Cˢⁱᵗₖ₁(Θᴷ,Sᵖ,Pˢⁱᵗₖ₁)

Return the first dissociation constant of silicic acid (H4SiO4) in seawater, given 
temperature in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pˢⁱᵗₖ₁`.

References: Yao and Millero (1995) cited by Millero (1995)
pH scale  : SWS (according to Dickson et al, 2007)
Note      : No pressure correction available
Note      : converted here from mol/kg-H2O to mol/kg-sw
"""
function Fˢⁱᵗₖ₁(Θᴷ,Sᵖ,Pˢⁱᵗₖ₁)
    ln_ksi1_p0 = ( Pˢⁱᵗₖ₁.a₀ +
                   Pˢⁱᵗₖ₁.a₁/Θᴷ +
                   Pˢⁱᵗₖ₁.a₂*log(Θᴷ) +
                  ( Pˢⁱᵗₖ₁.b₀ +
                    Pˢⁱᵗₖ₁.b₁/Θᴷ
                  )*sqrt(μₛ(Sᵖ)) +
                  ( Pˢⁱᵗₖ₁.b₂ +
                    Pˢⁱᵗₖ₁.b₃/Θᴷ
                  )*μₛ(Sᵖ) +
                  ( Pˢⁱᵗₖ₁.c₀ +
                    Pˢⁱᵗₖ₁.c₁/Θᴷ
                  )*μₛ(Sᵖ)*μₛ(Sᵖ) )

    # Pressure correction : currently none
    ln_ksi1_pp = 0.

    exp( ln_ksi1_p0 + ln_ksi1_pp ) * H₂Oˢʷ(Sᵖ)
end

Base.@kwdef struct Pᴴ²ˢₖ₁{FT}
    a₀ :: FT =     225.838
    a₁ :: FT = - 13275.3
    a₂ :: FT = -    34.6435
    a₃ :: FT =       0.3449
    a₄ :: FT = -     0.0274
    v₀ :: FT = -    14.80
    v₁ :: FT =       0.0020
    v₂ :: FT =  -    0.400E-03
    k₀ :: FT =       2.89e-3
    k₁ :: FT =       0.054e-3
end

@inline
"""
Cᴴ²ˢₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴ²ˢₖ₁)

Returns the dissociation constant of hydrogen sulfide in sea-water, given temperature in K,
`Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᴴ²ˢₖ₁`.

References: Millero et al. (1988) (cited by Millero (1995)
            Millero (1995) for pressure correction
pH scale  : - SWS (according to Yao and Millero, 1995, p. 82: "refitted if necessary")
            - Total (according to Lewis and Wallace, 1998)
Note      : we stick to SWS here for the time being
Note      : the fits from Millero (1995) and Yao and Millero (1995)
            derive from Millero et al. (1988), with all the coefficients
            multiplied by -ln(10)
"""
function Fᴴ²ˢₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴ²ˢₖ₁)
    ln_kh2s_p0 = ( Pᴴ²ˢₖ₁.a₀ +
                   Pᴴ²ˢₖ₁.a₁/Θᴷ +
                   Pᴴ²ˢₖ₁.a₂*log(Θᴷ) +
                   Pᴴ²ˢₖ₁.a₃*Sᴾ⁰⁵(Sᵖ) +
                   Pᴴ²ˢₖ₁.a₄*Sᵖ)

    "Pressure correction for applied pressure /= 0"
    ln_kh2s_pp = (-( Pᴴ²ˢₖ₁.v₀ +
                     Pᴴ²ˢₖ₁.v₁*ΘᵒC(Θᴷ) +
                     Pᴴ²ˢₖ₁.v₂*ΘᵒC(Θᴷ)*ΘᵒC(Θᴷ)
                   )+( 
                      Pᴴ²ˢₖ₁.k₀ +
                      Pᴴ²ˢₖ₁.k₁*ΘᵒC(Θᴷ)
                     )*(Δpᵦₐᵣ/2.)
                     )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    exp( ln_kh2s_p0 + ln_kh2s_pp )
end

Base.@kwdef struct Pᴺᴴ⁴ₖ₁{FT}
    a₀ :: FT = -    0.25444
    a₁ :: FT = - 6285.33
    a₂ :: FT =      0.0001635
    b₀ :: FT =      0.46532
    b₁ :: FT = -  123.7184
    b₂ :: FT = -    0.01992
    b₃ :: FT =      3.17556
    v₀ :: FT = -   26.43
    v₁ :: FT =      0.0889
    v₂ :: FT = -    0.905E-03
    k₀ :: FT = -    5.03E-03
    k₁ :: FT =      0.0814E-03
end

@inline
"""
Cᴺᴴ⁴ₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴺᴴ⁴ₖ₁)

Return the dissociation constant of ammonium in sea-water [mol/kg-SW], given temperature
in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᴺᴴ⁴ₖ₁`.

References: Yao and Millero (1995)
            Millero (1995) for pressure correction
pH scale  : SWS
"""
function Fᴺᴴ⁴ₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴺᴴ⁴ₖ₁)
    ln_knh4_p0 = ( Pᴺᴴ⁴ₖ₁.a₀ +
                     Pᴺᴴ⁴ₖ₁.a₁/Θᴷ +
                        Pᴺᴴ⁴ₖ₁.a₂*Θᴷ +
                     ( Pᴺᴴ⁴ₖ₁.b₀ +
                          Pᴺᴴ⁴ₖ₁.b₁/Θᴷ
                     )*Sᴾ⁰⁵(Sᵖ) +
                     ( Pᴺᴴ⁴ₖ₁.b₂ +
                          Pᴺᴴ⁴ₖ₁.b₃/Θᴷ
                     )*Sᵖ)

    
    "Pressure correction for applied pressure /= 0"
    ln_knh4_pp = (-( Pᴺᴴ⁴ₖ₁.v₀ +
                     Pᴺᴴ⁴ₖ₁.v₁*ΘᵒC(Θᴷ) +
                     Pᴺᴴ⁴ₖ₁.v₂*ΘᵒC(Θᴷ)*ΘᵒC(Θᴷ)
                   )+( 
                      Pᴺᴴ⁴ₖ₁.k₀ +
                      Pᴺᴴ⁴ₖ₁.k₁*ΘᵒC(Θᴷ)
                     )*(Δpᵦₐᵣ/2.)
                     )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    exp( ln_knh4_p0 + ln_knh4_pp )
end

Base.@kwdef struct Pᴴᶠᵦ₁{FT}
    a₀ :: FT =     12.641
    a₁ :: FT = - 1590.2
    a₂ :: FT = -    1.525
    v₀ :: FT = -    9.78
    v₁ :: FT = -    0.0090
    v₂ :: FT = -    0.942E-03
    k₀ :: FT = -    3.91e-3
    k₁ :: FT =      0.054e-3
end

@inline
"""
Cᴴᶠᵦ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴᶠᵦ₁)

Return the association constant of HF in (mol/kg-SW)^-1, , given temperature
in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᴴᶠᵦ₁`.

HF <-> H⁺ + F⁻ 

References: Dickson and Riley (1979)
            Millero (1995) for pressure correction
pH scale  : free
Note      : converted here from mol/kg-H2O to mol/kg-SW
"""
function Fᴴᶠᵦ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴᶠᵦ₁)
    ln_bhf_p0 = ( Pᴴᶠᵦ₁.a₀ +
                  Pᴴᶠᵦ₁.a₁/Θᴷ +
                  Pᴴᶠᵦ₁.a₂*sqrt(μₛ(Sᵖ)))

    "Pressure correction for applied pressure /= 0"
    ln_khf_pp = (-( Pᴴᶠᵦ₁.v₀ +
                     Pᴴᶠᵦ₁.v₁*ΘᵒC(Θᴷ) +
                     Pᴴᶠᵦ₁.v₂*ΘᵒC(Θᴷ)*ΘᵒC(Θᴷ)
                   )+( 
                      Pᴴᶠᵦ₁.k₀ +
                      Pᴴᶠᵦ₁.k₁*ΘᵒC(Θᴷ)
                     )*(Δpᵦₐᵣ/2.)
                     )*(Δpᵦₐᵣ/Rₜ(Θᴷ))
    """                
    Final \beta_HF value
    --------------------
     notice that  ln(k_HF(P)) = ln(k_HF(0)) + zln_khf_pp
            <=>  -ln(\beta_HF(P)) = -ln(\beta_HF(0)) + zln_khf_pp
            <=>   ln(\beta_HF(P)) =  ln(\beta_HF(0)) - zln_khf_pp
    """
    exp(ln_bhf_p0 - ln_khf_pp ) / H₂Oˢʷ(Sᵖ)
end

Base.@kwdef struct Pᴴᶠₖ₁{FT}
    a₀ :: FT = -  9.68
    a₁ :: FT =  874.0
    a₂ :: FT =    0.111
    v₀ :: FT = -  9.78
    v₁ :: FT = -  0.0090
    v₂ :: FT = -  0.942E-3
    k₀ :: FT = -  3.91e-3
    k₁ :: FT =    0.054e-3
end

@inline """
Cᴴᶠₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴᶠₖ₁)

Return the dissociation constant for hydrogen fluoride in mol/kg-SW, given temperature
in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᴴᶠₖ₁`.

HF <-> H⁺ + F⁻ 

References: Perez and Fraga (1987)
            Millero (1995) for pressure correction
pH scale  : Total (according to Handbook, 2007
"""
function Fᴴᶠₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴᶠₖ₁)
    ln_khf_p0 = ( Pᴴᶠₖ₁.a₀ +
                  Pᴴᶠₖ₁.a₁/Θᴷ +
                  Pᴴᶠₖ₁.a₂*Sᴾ⁰⁵(Sᵖ))

    "Pressure correction for applied pressure /= 0"
    ln_khf_pp = (-( Pᴴᶠₖ₁.v₀ +
                     Pᴴᶠₖ₁.v₁*ΘᵒC(Θᴷ) +
                     Pᴴᶠₖ₁.v₂*ΘᵒC(Θᴷ)*ΘᵒC(Θᴷ)
                   )+( 
                      Pᴴᶠₖ₁.k₀ +
                      Pᴴᶠₖ₁.k₁*ΘᵒC(Θᴷ)
                     )*(Δpᵦₐᵣ/2.)
                     )*(Δpᵦₐᵣ/Rₜ(Θᴷ))
        
    exp( ln_khf_p0 + ln_khf_pp )
end

Base.@kwdef struct Pᴴˢᴼ⁴ₖ₁{FT}
    a₀ :: FT =     141.328
    a₁ :: FT = -  4276.1
    a₂ :: FT = -    23.093
    b₀ :: FT =     324.57
    b₁ :: FT = - 13856.
    b₂ :: FT = -    47.986
    c₀ :: FT = -   771.54
    c₁ :: FT =   35474.
    c₂ :: FT =     114.723
    d₀ :: FT = -  2698.
    d₁ :: FT =    1776.
    v₀ :: FT = -    18.03
    v₁ :: FT =       0.0466
    v₂ :: FT =       0.316E-03
    k₀ :: FT = -     4.53e-3
    k₁ :: FT =       0.0900e-3
end

@inline """
Cᴴˢᴼ⁴ₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴˢᴼ⁴ₖ₁)

Function returns the dissociation constant of hydrogen sulfate (bisulfate) , given temperature
in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᴴˢᴼ⁴ₖ₁`.

References: Dickson (1990) -- also Handbook (2007)
            Millero (1995) for pressure correction
pH scale  : free
Note      : converted here from mol/kg-H2O to mol/kg-SW
"""
function Fᴴˢᴼ⁴ₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴˢᴼ⁴ₖ₁)
    ln_khso4_p0 = ( Pᴴˢᴼ⁴ₖ₁.a₀ +
                    Pᴴˢᴼ⁴ₖ₁.a₁/Θᴷ +
                    Pᴴˢᴼ⁴ₖ₁.a₂*log(Θᴷ) +
                   (Pᴴˢᴼ⁴ₖ₁.b₀ +
                    Pᴴˢᴼ⁴ₖ₁.b₁/Θᴷ +
                    Pᴴˢᴼ⁴ₖ₁.b₂*log(Θᴷ)
                   )*sqrt(μₛ(Sᵖ)) +
                   (Pᴴˢᴼ⁴ₖ₁.c₀ +
                    Pᴴˢᴼ⁴ₖ₁.c₁/Θᴷ +
                    Pᴴˢᴼ⁴ₖ₁.c₂*log(Θᴷ)
                   )*μₛ(Sᵖ) +
                   (Pᴴˢᴼ⁴ₖ₁.d₀/Θᴷ)*sqrt(μₛ(Sᵖ))*μₛ(Sᵖ)+
                   (Pᴴˢᴼ⁴ₖ₁.d₁/Θᴷ)*μₛ(Sᵖ)*μₛ(Sᵖ))

    "Pressure correction for applied pressure /= 0"
    ln_khso4_pp = (-( Pᴴˢᴼ⁴ₖ₁.v₀ +
                      Pᴴˢᴼ⁴ₖ₁.v₁*ΘᵒC(Θᴷ) +
                      Pᴴˢᴼ⁴ₖ₁.v₂*ΘᵒC(Θᴷ)*ΘᵒC(Θᴷ)
                    )+( 
                       Pᴴˢᴼ⁴ₖ₁.k₀ +
                       Pᴴˢᴼ⁴ₖ₁.k₁*ΘᵒC(Θᴷ)
                      )*(Δpᵦₐᵣ/2.)
                      )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    exp( ln_khso4_p0 + ln_khso4_pp ) * H₂Oˢʷ(Sᵖ)
end

Base.@kwdef struct Pᶜᵃˡᶜⁱᵗᵉₛₚ{FT}
    a₀ :: FT = - 171.9065
    a₁ :: FT = -   0.077993
    a₂ :: FT =  2839.319
    a₃ :: FT =    71.595
    b₀ :: FT = -   0.77712
    b₁ :: FT =     0.0028426
    b₂ :: FT =   178.34
    c₀ :: FT = -   0.07711
    d₀ :: FT =     0.0041249
    v₀ :: FT = -  48.76
    v₁ :: FT =     0.5304
    k₀ :: FT = -  11.76e-3
    k₁ :: FT =     0.3692e-3
end

@inline
"""
    Cᶜᵃˡᶜⁱᵗᵉₛₚ(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᶜᵃˡᶜⁱᵗᵉₛₚ)

Return the stoechiometric solubility product of calcite, `Ω`, in seawater, given temperature
in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᶜᵃˡᶜⁱᵗᵉₛₚ`

References: Mucci (1983)
            Millero (1995) for pressure correction
pH scale  : N/A
Units     : (mol/kg-SW)^2
"""
function Fᶜᵃˡᶜⁱᵗᵉₛₚ(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᶜᵃˡᶜⁱᵗᵉₛₚ)
    log10_kcalcite_p0 = ( Pᶜᵃˡᶜⁱᵗᵉₛₚ.a₀ +
                          Pᶜᵃˡᶜⁱᵗᵉₛₚ.a₁*Θᴷ +
                          Pᶜᵃˡᶜⁱᵗᵉₛₚ.a₂/Θᴷ +
                          Pᶜᵃˡᶜⁱᵗᵉₛₚ.a₃*log10(Θᴷ) +
                         (Pᶜᵃˡᶜⁱᵗᵉₛₚ.b₀ +
                          Pᶜᵃˡᶜⁱᵗᵉₛₚ.b₁*Θᴷ +
                          Pᶜᵃˡᶜⁱᵗᵉₛₚ.b₂/Θᴷ
                         )*Sᴾ⁰⁵(Sᵖ) +
                         (Pᶜᵃˡᶜⁱᵗᵉₛₚ.c₀)*Sᵖ +
                         (Pᶜᵃˡᶜⁱᵗᵉₛₚ.d₀)*Sᵖ*Sᴾ⁰⁵(Sᵖ) )

    "Pressure correction for applied pressure /= 0"
    ln_kcalcite_pp = (-( Pᶜᵃˡᶜⁱᵗᵉₛₚ.v₀ +
                            Pᶜᵃˡᶜⁱᵗᵉₛₚ.v₁*ΘᵒC(Θᴷ)
                          )+( 
                             Pᶜᵃˡᶜⁱᵗᵉₛₚ.k₀ +
                             Pᶜᵃˡᶜⁱᵗᵉₛₚ.k₁*ΘᵒC(Θᴷ)
                            )*(Δpᵦₐᵣ/2.)
                            )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

 10^(log10_kcalcite_p0) * exp(ln_kcalcite_pp)
end

Base.@kwdef struct Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ{FT}
    a₀ :: FT = - 171.945
    a₁ :: FT = -   0.077993
    a₂ :: FT =  2903.293
    a₃ :: FT =    71.595
    b₀ :: FT = -   0.068393
    b₁ :: FT =     0.0017276
    b₂ :: FT =    88.135
    c₀ :: FT = -   0.10018
    d₀ :: FT =     0.0059415
    v₀ :: FT = -  48.76
    v₁ :: FT =     0.5304
    v₂ :: FT =     2.8
    k₀ :: FT = -  11.76e-3
    k₁ :: FT =     0.3692e-3
end

@inline
"""
Cᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ)

Function returns stoechiometric solubility product, `Ω`, of aragonite in seawater, given temperature
in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ`.

References: Mucci (1983)
            Millero (1979) for pressure correction
pH scale  : N/A
Units     : (mol/kg-SW)^2
"""
function Fᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ)
    log10_karagonite_p0 = ( Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ.a₀ +
                            Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ.a₁*Θᴷ +
                            Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ.a₂/Θᴷ +
                            Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ.a₃*log10(Θᴷ) +
                           (Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ.b₀ +
                            Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ.b₁*Θᴷ +
                            Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ.b₂/Θᴷ
                           )*Sᴾ⁰⁵(Sᵖ) +
                           (Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ.c₀)*Sᵖ +
                           (Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ.d₀)*Sᵖ*Sᴾ⁰⁵(Sᵖ) )

    "Pressure correction for applied pressure /= 0"
    ln_karagonite_pp = (-( Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ.v₀ +
                           Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ.v₁*ΘᵒC(Θᴷ) +
                           Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ.v₂
                          )+( 
                             Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ.k₀ +
                             Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ.k₁*ΘᵒC(Θᴷ)
                            )*(Δpᵦₐᵣ/2.)
                            )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    10^(log10_karagonite_p0) * exp(ln_karagonite_pp)
end

@inline
"""
H⁺ₛoverH⁺ₜ(Θᴷ,Sᵖ,Δpᵦₐᵣ)

Returns the ratio H_SWS/H_Tot as a function of salinity, `Sᵖ`.
Reference:  Munhoven
pH scale:   all
"""
function H⁺ₛoverH⁺ₜ(Θᴷ,Sᵖ,Δpᵦₐᵣ)
    (1. +  
            (Fᵀᴼᵀ(Sᵖ)*Fᴴᶠᵦ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴᶠᵦ₁))
           /(1. + SO₄ᵀᴼᵀ(Sᵖ)/Fᴴˢᴼ⁴ₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴˢᴼ⁴ₖ₁))
           )
end

@inline
"""
H⁺ₜoverH⁺₃(Θᴷ,Sᵖ,Δpᵦₐᵣ)

Function returns the ratio H_Tot/H_free as a function of salinity, `Sᵖ`.
Reference:  Munhoven
pH scale:   N/A
"""
function H⁺ₜoverH⁺₃(Θᴷ,Sᵖ,Δpᵦₐᵣ)
    1. + SO₄ᵀᴼᵀ(Sᵖ)/Fᴴˢᴼ⁴ₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴˢᴼ⁴ₖ₁)
end

@inline
"""
H⁺ₛoverH⁺₃(Θᴷ,Sᵖ,Δpᵦₐᵣ)

Function returns the ratio H_SWS/H_free as a function
of salinity, `Sᵖ`.
"""
function H⁺ₛoverH⁺₃(Θᴷ,Sᵖ,Δpᵦₐᵣ)
    ( 1. + 
            (Fᵀᴼᵀ(Sᵖ)*Fᴴᶠᵦ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴᶠᵦ₁)) +
            (SO₄ᵀᴼᵀ(Sᵖ)/Fᴴˢᴼ⁴ₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴˢᴼ⁴ₖ₁))
         )
end
"""
   Copyright  2013, 2014, 2020, 2021 Guy Munhoven
   This file is part of SolveSAPHE v. 2
   SolveSAPHE is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   SolveSAPHE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.
   You should have received a copy of the GNU Lesser General Public License
   along with SolveSAPHE.  If not, see <http://www.gnu.org/licenses/>.
"""
end # module DissociationConstants
#------------------------------------------------------------