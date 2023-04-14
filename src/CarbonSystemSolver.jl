module CarbonSystemSolver
#
#    Copyright  2013, 2014, 2020, 2021 Guy Munhoven
#
#    This file is part of SolveSAPHE v. 2
#
#    SolveSAPHE is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SolveSAPHE is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with SolveSAPHE.  If not, see <http://www.gnu.org/licenses/>.
#

# ------------------------------------------------------------
module DissociationConstants
export DissociationCoefficients

struct DissociationCoefficients
#    {FT}
    """
    DissociationCoefficients(FT=Float64)
    Function to calculate dissociation coefficients of the carbon system.
    """
#    FT    :: Float64
#    "Temperature in degrees Celsius"
#    Θ     :: FT 
#    "Absolute Salinity in g/kg"
#    Sᴬ    :: FT
#    "Applied Pressure in bars"
#    Δpᵦₐᵣ :: FT
    "Temperature in degrees Celsius"
    Θ     :: Float64
    "Absolute Salinity in g/kg"
    Sᴬ    :: Float64
    "Applied Pressure in bars"
    Δpᵦₐᵣ :: Float64
    function DissociationCoefficients(
#        FT          = Float64,
#        Θ     :: FT = 25.0,
#        Sᴬ    :: FT = 35.0,
#        Δpᵦₐᵣ :: FT = 0.0,
#        FT          = Float64,
        Θᶜ    :: Float64 = 25.0,
        Sᴬ    :: Float64 = 35.0,
        Δpᵦₐᵣ :: Float64 = 0.0,
        )

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

@inline Θᴷ₁₀₀(Θᴷ) = Θᴷ/100
@inline Sᴾ⁰⁵(Sᵖ)  = sqrt(Sᵖ)
@inline ΘᵒC(Θᴷ)   = Θᴷ - Θᴷ_0ᵒC
@inline ΔSᵖ(Sᵖ)   = Sᵖ - 34.8
@inline Rₜ(Θᴷ)     = gasconst_bar_cm3_o_mol_k * Θᴷ

@inline """
    H₂Oˢʷ(Sᵖ)

    # Function returns the mass of pure water in one kg of seawater
    # of salinity s
    # References: "libthdyct" -- derived by Munhoven (1997) from data by Millero (1982)
    #             "Handbook (2007)" -- Handbook (2007)
    # pH scale:   N/A
"""
H₂Oˢʷ(Sᵖ) =  1. - 0.001005*Sᵖ 
#H₂Oˢʷ(Sᵖ) = 1. - 0.0010049*Sᵖ # libthdyct

@inline """
    μₛ(Sᵖ)

    # Function calculates ionic strength in mol/kg-SW, for given salinity.
    # References: "libthdyct" -- derived by Munhoven (1997) from data by Millero (1982)
    #             "Handbook (2007)" -- Handbook (2007)
    # pH scale:   N/A
"""
μₛ(Sᵖ)    = (0.019924 * Sᵖ) / H₂Oˢʷ(Sᵖ) # Handbook (2007)
# μₛ(Sᵖ)    = (0.019920 * Sᵖ) / H₂Oˢʷ(Sᵖ)# libthdyct

@inline """
    Bᵀᴼᵀ(Sᵖ)

    # Total borate concentration in mol/kg-SW given the salinity of a sample
    # References: Uppström (1974), cited by  Dickson et al. (2007, chapter 5, p 10)
    #             Millero (1982) cited in Millero (1995)
"""
Bᵀᴼᵀ(Sᵖ)  = 0.000416*(Sᵖ/35.)
#Bᵀᴼᵀ(Sᵖ)  = 0.000232 * (Sᵖ/1.80655)/10.811

@inline """
    Caᵀᴼᵀ(Sᵖ)

    # Total calcium concentration in mol/kg-SW given the salinity of a sample
    # References: Culkin and Cox (1966), Culkin (1967), Riley and Tongudai (1967)    
"""
Caᵀᴼᵀ(Sᵖ) = (0.02127/40.078) * (Sᵖ/1.80655)
# Caᵀᴼᵀ(Sᵖ) = 0.010282*(Sᵖ/35.)

@inline """
    Fᵀᴼᵀ(Sᵖ)

    # Total fluoride concentration in mol/kg-SW given the salinity of a sample
    # References: Culkin (1965) (???)
"""
Fᵀᴼᵀ(Sᵖ)  = 0.000068*(Sᵖ/35.)

@inline """
    SO₄ᵀᴼᵀ(Sᵖ)

    # Total sulfate concentration in mol/kg-SW given the salinity of a sample
    # References: Morris, A.W. and Riley, J.P. (1966) quoted in Handbook (2007)
"""
SO₄ᵀᴼᵀ(Sᵖ) = (0.1400/96.062)*(Sᵖ/1.80655)
#SO₄ᵀᴼᵀ(Sᵖ) = 0.028234*(Sᵖ/35.)

@inline """
    NH₄ᵀᴼᵀ(Sᵖ)

    # Total ammonium concentration in mol/kg-SW given the salinity of a sample
    # References: Yamamoto (1995)
"""
NH₄ᵀᴼᵀ(Sᵖ) = (0.0001/17.031) * (Sᵖ/1.80655)

@inline """
    PO₄ᵀᴼᵀ(Sᵖ)

    # Total phosphate concentration in mol/kg-SW given the salinity of a sample
    # References: Millero (1995)
"""
PO₄ᵀᴼᵀ(Sᵖ) = (0.0001/94.971) * (Sᵖ/1.80655)

@inline """
    SiO₃ᵀᴼᵀ(Sᵖ)

    # Total silicate concentration in mol/kg-SW given the salinity of a sample
    # References: Millero (1995)
"""
SiO₃ᵀᴼᵀ(Sᵖ) = (0.0001/60.084) * (Sᵖ/1.80655)

@inline """
    H₂Sᵀᴼᵀ(Sᵖ)

    # Total hydrogen sulfide concentration in mol/kg-SW given the salinity of a sample
    # References: Dickson (1990)
"""
H₂Sᵀᴼᵀ(Sᵖ) = (0.0001/34.082) * (Sᵖ/1.80655)

##=======================================================================
## AK_CARB_0_WEIS74 <- function (Θᴷ, s)
##=======================================================================
const Pᵈⁱᶜₖ₀ = (
    a₀ = -60.2409,
    a₁ =  93.4517,
    a₂ =  23.3585,
    b₀ =   0.023517,
    b₁ = - 0.023656,
    b₂ =   0.0047036,
)

@inline """
    Cᵈⁱᶜₖ₀(Θᴷ,Sᵖ,Pᵈⁱᶜₖ₀)

    # Function calculates K0 in (mol/kg-SW)/atmosphere
    # References: Weiss (1979) [(mol/kg-SW)/atm]
    # pH scale  : N/A
    # Note      : currently no pressure correction
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
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

##=======================================================================
##AK_CARB_1_ROYE93 <- function (Θᴷ, s, p_bar)
##=======================================================================
const Pᵈⁱᶜₖ₁ᵣ₉₃ = (
    a₀ =     2.83655,
    a₁ = -2307.1266,
    a₂ = -   1.5529413,
    b₀ = -   0.20760841,
    b₁ = -   4.0484,
    b₂ =     0.08468345,
    b₃ = -   0.00654208,
    v₀ = -  25.5,
    v₁ = -   0.151,
    v₂ =     0.1271,
    k₀ = -   3.08e-3,
    k₁ = -   0.578e-3,
    k₂ =     0.0877e-3,
)

@inline """
    Cᵈⁱᶜₖ₁ᵣ₉₃(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᵈⁱᶜₖ₁ᵣ₉₃)

    # Function calculates first dissociation constant of carbonic acid
    # in mol/kg-SW on the Total pH-scale.
    # References: Roy et al. (1993) -- also Handbook (1994)
    #             Millero (1979) pressure correction
    # pH scale  : Total
    # Note      : converted here from mol/kg-H2O to mol/kg-SW
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
    #     Δpᵦₐᵣ  : applied pressure in bar
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

    # Pressure correction for applied pressure p_bar /= 0
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
 
##=======================================================================
##AK_CARB_2_ROYE93 <- function (Θᴷ, s, p_bar)
##=======================================================================
const Pᵈⁱᶜₖ₂ᵣ₉₃ = (
    a₀ = -   9.226508,
    a₁ = -3351.6106,
    a₂ = -   0.2005743,
    b₀ = -   0.106901773,
    b₁ = -  23.9722,
    b₂ =     0.1130822,
    b₃ = -   0.00846934,
    v₀ = -  15.82,
    v₁ =     0.321,
    v₂ = -   0.0219,
    k₀ =     1.13e-3,
    k₁ = -   0.314e-3,
    k₂ = -   0.1475e-3,
)

@inline """
    Cᵈⁱᶜₖ₂ᵣ₉₃(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᵈⁱᶜₖ₂ᵣ₉₃)

    # Function calculates second dissociation constant K1
    # in mol/kg-SW on the Total pH-scale.
    # References: Roy et al. (1993) -- also Handbook (1994)
    #             Millero (1979) pressure correction
    # pH scale  : Total
    # Note      : converted here from mol/kg-H2O to mol/kg-SW
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
    #     Δpᵦₐᵣ  : applied pressure in bar
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

    # Pressure correction for applied pressure p_bar /= 0
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

##=======================================================================
##AK_CARB_1_MILL95 <- function (Θᴷ, s, p_bar)
##=======================================================================
const Pᵈⁱᶜₖ₁ₘ₉₅ = (
    a₀ =     2.18867,
    a₁ = -2275.0360,
    a₂ = -   1.468591,
    b₀ = -   0.138681,
    b₁ = -   9.33291,
    b₂ =     0.0726483,
    b₃ = -   0.00574938,
    v₀ = -  25.5,
    v₁ = -   0.151,
    v₂ =     0.1271,
    k₀ = -   3.08e-3,
    k₁ = -   0.578e-3,
    k₂ =     0.0877e-3,
)

@inline """
    Cᵈⁱᶜₖ₁ₘ₉₅(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᵈⁱᶜₖ₁ₘ₉₅)

    # Function calculates first dissociation constant of carbonic acid
    # in mol/kg-SW on the SWS pH-scale.
    # References: Millero (1995, eq 50 -- ln K1(COM))
    #             Millero (1982) pressure correction
    # pH scale:   SWS
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
    #     Δpᵦₐᵣ  : applied pressure in bar
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


    # Pressure correction for applied pressure p_bar /= 0
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

##=======================================================================
##AK_CARB_2_MILL95 <- function (Θᴷ, s, p_bar)
##=======================================================================
const Pᵈⁱᶜₖ₂ₘ₉₅ = (
    a₀ = -   0.84226,
    a₁ = -3741.1288,
    a₂ = -   1.437139,
    b₀ = -   0.128417,
    b₁ = -  24.41239,
    b₂ =     0.1195308,
    b₃ = -   0.00912840,
    v₀ = -  15.82,
    v₁ =     0.321,
    v₂ = -   0.0219,
    k₀ =     1.13e-3,
    k₁ = -   0.314e-3,
    k₂ = -   0.1475e-3,
)

@inline """
    Cᵈⁱᶜₖ₂ₘ₉₅(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᵈⁱᶜₖ₂ₘ₉₅)

    # Function calculates second dissociation constant K1
    # in mol/kg-SW on the SWS pH-scale.
    # References: Millero (1995, eq 51 -- ln K2(COM))
    #             Millero (1979) pressure correction
    # pH scale:   SWS
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
    #     Δpᵦₐᵣ  : applied pressure in bar
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

    # Pressure correction for applied pressure p_bar /= 0
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

##=======================================================================
##AK_CARB_1_LUEK00 <- function (Θᴷ, s, p_bar)
##=======================================================================
const Pᵈⁱᶜₖ₁ₗ₀₀ = (
    a₀ =    61.2172,
    a₁ = -3633.86,
    a₂ = -   9.67770,
    b₀ =     0.011555,
    b₁ = -   0.0001152,
    v₀ = -  25.5,
    v₁ = -   0.151,
    v₂ =     0.1271,
    k₀ = -   3.08e-3,
    k₁ = -   0.578e-3,
    k₂ =     0.0877e-3,
)

@inline """
    Cᵈⁱᶜₖ₁ₗ₀₀(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᵈⁱᶜₖ₁ₗ₀₀)

    # Function calculates first dissociation constant of carbonic acid
    # in mol/kg-SW on the Total pH-scale.
    # References: Luecker et al. (2000) -- also Handbook (2007)
    #             Millero (1979) pressure correction
    # pH scale:   Total
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
    #     Δpᵦₐᵣ  : applied pressure in bar
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

    # Pressure correction for applied pressure p_bar /= 0
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

##=======================================================================
##AK_CARB_2_LUEK00 <- function (Θᴷ, s, p_bar)
##=======================================================================
const Pᵈⁱᶜₖ₂ₗ₀₀ = (
    a₀ = -  25.9290,
    a₁ = - 471.78,
    a₂ =     3.16967,
    b₀ =     0.01781,
    b₁ = -   0.0001122,
    v₀ = -  15.82,
    v₁ =     0.321,
    v₂ = -   0.0219,
    k₀ =     1.13e-3,
    k₁ = -   0.314e-3,
    k₂ = -   0.1475e-3,
)

@inline """
    Cᵈⁱᶜₖ₂ₗ₀₀(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᵈⁱᶜₖ₂ₗ₀₀)

    # Function calculates second dissociation constant K1
    # in mol/kg-SW on the Total pH-scale.
    # References: Luecker et al. (2000) -- also Handbook (2007)
    #             Millero (1979) pressure correction
    # pH scale:   Total
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
    #     Δpᵦₐᵣ  : applied pressure in bar
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

    # Pressure correction for applied pressure p_bar /= 0
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

##=======================================================================
##AK_BORA_DICK90 <- function (Θᴷ, s, p_bar)
##=======================================================================
const Pᴮₖ₁ = (
    a₀ = - 8966.90,
    a₁ = - 2890.53,
    a₂ = -   77.942,
    a₃ =      1.728,
    a₄ = -    0.0996,
    b₀ =    148.0248,
    b₁ =    137.1942,
    b₂ =      1.62142,
    c₀ = -   24.4344,
    c₁ = -   25.085,
    c₂ = -    0.2474,
    d₀ =      0.053105,
    v₀ = -   29.48,
    v₁ =      0.295,
    v₂ =      0.1622,
    v₃ = -    0.002608,
    k₀ = -    2.84e-3,
    k₁ =      0.354e-3,
)

@inline """
    Cᵇₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴮₖ₁)

    # Function calculates boric acid dissociation constant KB
    # in mol/kg-SW on the total pH-scale.
    # References: Dickson (1990, eq. 23) -- also Handbook (2007, eq. 37)
    #             Millero (1979) pressure correction
    # pH scale  : total
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
    #     Δpᵦₐᵣ  : applied pressure in bar
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

    # Pressure correction for applied pressure p_bar /= 0
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

##=======================================================================
##AK_W_MILL95 <- function (Θᴷ, s, p_bar)
##=======================================================================
const Pᴴ²ᴼₖ₁ = (
    a₀ =    148.9802,
    a₁ = -13847.26,
    a₂ = -   23.6521,
    b₀ = -    5.977,
    b₁ =    118.67,
    b₂ =      1.0495,
    c₀ = -    0.01615,
    v₀ =  -  20.02,
    v₁ =      0.1119,
    v₂ =  -   0.1409E-02,
    k₀ = -    5.13e-3,
    k₁ =      0.0794e-3
)

@inline """
    Cᴴ²ᴼₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴ²ᴼₖ₁)

    # Function calculates water dissociation constant Kw in (mol/kg-SW)^2
    # References: Millero (1995) for value at p_bar = 0
    #             Millero (pers. comm. 1996) for pressure correction
    # pH scale  : SWS
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
    #     Δpᵦₐᵣ  : applied pressure in bar
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

    # Pressure correction for applied pressure p_bar /= 0
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

##=======================================================================
##AK_PHOS_1_MILL95 <- function (Θᴷ, s, p_bar)
##=======================================================================
const Pᴾᴼ⁴ₖ₁ = (
    a₀ =   115.54,
    a₁ = -4576.752,
    a₂ = -  18.453,
    b₀ =     0.69171,
    b₁ = - 106.736,
    b₂ = -   0.01844,
    b₃ = -   0.65643,
    v₀ = -  14.51,
    v₁ =     0.1211,
    v₂ = -   0.321E-03,
    k₀ = -   2.67e-3,
    k₁ =     0.0427e-3
)
@inline """
    Cᴾᴼ⁴ₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴾᴼ⁴ₖ₁)

    # Function returns the first dissociation constant
    # of phosphoric acid (H3PO4) in seawater
    # References: Yao and Millero (1995)
    #             Millero (1995) for pressure correction
    # pH scale  : SWS
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
    #     Δpᵦₐᵣ  : applied pressure in bar
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

    # Pressure correction for applied pressure p_bar /= 0
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

##=======================================================================
##AK_PHOS_2_MILL95 <- function (Θᴷ, s, p_bar)
##=======================================================================
const Pᴾᴼ⁴ₖ₂ = (
    a₀ =    172.1033,
    a₁ = - 8814.715,
    a₂ = -   27.927,
    b₀ =      1.3566,
    b₁ = -  160.340,
    b₂ = -    0.05778,
    b₃ =      0.37335,
    v₀ = -   23.12,
    v₁ =      0.1758,
    v₂ = -    0.002647,
    k₀ = -    5.15e-3,
    k₁ =      0.09e-3
)

@inline """
    Cᴾᴼ⁴ₖ₂(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴾᴼ⁴ₖ₂)

    # Function returns the second dissociation constant
    # of phosphoric acid (H3PO4) in seawater
    # References: Yao and Millero (1995)
    #             Millero (1995) for pressure correction
    # pH scale  : SWS
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
    #     Δpᵦₐᵣ  : applied pressure in bar
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

    # Pressure correction for applied pressure p_bar /= 0
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

##=======================================================================
##AK_PHOS_3_MILL95 <- function (Θᴷ, s, p_bar)
##=======================================================================
const Pᴾᴼ⁴ₖ₃ = (
    a₀ = -   18.126,
    a₁ = - 3070.75,
    a₂ =      2.81197,
    a₃ =     17.27039,
    a₄ = -    0.09984,
    a₅ = -   44.99486,
    v₀ = -   26.57,
    v₁ =      0.2020,
    v₂ = -    3.042e-3,
    k₀ = -    4.08e-3,
    k₁ =      0.0714e-3
)

@inline """
    Cᴾᴼ⁴ₖ₃(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴾᴼ⁴ₖ₃)

    # Function returns the third dissociation constant
    # of phosphoric acid (H3PO4) in seawater
    # References: Yao and Millero (1995)
    #             Millero (1995) for pressure correction
    # pH scale  : SWS
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
    #     Δpᵦₐᵣ  : applied pressure in bar
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

    # Pressure correction for applied pressure p_bar /= 0
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

##=======================================================================
##AK_SILI_1_MILL95 <- function (Θᴷ, s)
##=======================================================================
const Pˢⁱᵗₖ₁ = (
    a₀ =    117.40,
    a₁ = - 8904.2,
    a₂ = -   19.334,
    b₀ =      3.5913,
    b₁ = -  458.79,
    b₂ = -    1.5998,
    b₃ =    188.74,
    c₀ =      0.07871,
    c₁ = -   12.1652
)

@inline """
    Cˢⁱᵗₖ₁(Θᴷ,Sᵖ,Pˢⁱᵗₖ₁)

    # Function returns the first dissociation constant
    # of silicic acid (H4SiO4) in seawater
    # References: Yao and Millero (1995) cited by Millero (1995)
    # pH scale  : SWS (according to Dickson et al, 2007)
    # Note      : No pressure correction available
    # Note      : converted here from mol/kg-H2O to mol/kg-sw
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
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

##=======================================================================
##AK_H2S_1_MILL95 <- function (Θᴷ, s, p_bar)
##=======================================================================
const Pᴴ²ˢₖ₁ = (
    a₀ =     225.838,
    a₁ = - 13275.3,
    a₂ = -    34.6435,
    a₃ =       0.3449,
    a₄ = -     0.0274,
    v₀ = -    14.80,
    v₁ =       0.0020,
    v₂ =  -    0.400E-03,
    k₀ =       2.89e-3,
    k₁ =       0.054e-3
)

@inline """
    Cᴴ²ˢₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴ²ˢₖ₁)

    # Function returns the dissociation constant 
    # of hydrogen sulfide in sea-water
    # References: Millero et al. (1988) (cited by Millero (1995)
    #             Millero (1995) for pressure correction
    # pH scale  : - SWS (according to Yao and Millero, 1995, p. 82: "refitted if necessary")
    #             - Total (according to Lewis and Wallace, 1998)
    # Note      : we stick to SWS here for the time being
    # Note      : the fits from Millero (1995) and Yao and Millero (1995)
    #             derive from Millero et al. (1988), with all the coefficients
    #             multiplied by -ln(10)
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
    #     Δpᵦₐᵣ  : applied pressure in bar
"""
function Fᴴ²ˢₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴ²ˢₖ₁)
    ln_kh2s_p0 = ( Pᴴ²ˢₖ₁.a₀ +
                   Pᴴ²ˢₖ₁.a₁/Θᴷ +
                   Pᴴ²ˢₖ₁.a₂*log(Θᴷ) +
                   Pᴴ²ˢₖ₁.a₃*Sᴾ⁰⁵(Sᵖ) +
                   Pᴴ²ˢₖ₁.a₄*Sᵖ)

    # Pressure correction for applied pressure p_bar /= 0
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

##=======================================================================
##AK_AMMO_1_YAMI95 <- function (Θᴷ, s, p_bar)
##=======================================================================
const Pᴺᴴ⁴ₖ₁ = (
    a₀ = -    0.25444,
    a₁ = - 6285.33,
    a₂ =      0.0001635,
    b₀ =      0.46532,
    b₁ = -  123.7184,
    b₂ = -    0.01992,
    b₃ =      3.17556,
    v₀ = -   26.43,
    v₁ =      0.0889,
    v₂ = -    0.905E-03,
    k₀ = -    5.03E-03,
    k₁ =      0.0814E-03
)

@inline """
    Cᴺᴴ⁴ₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴺᴴ⁴ₖ₁)

    # Function returns the dissociation constant
    # of ammonium in sea-water [mol/kg-SW]
    # References: Yao and Millero (1995)
    #             Millero (1995) for pressure correction
    # pH scale  : SWS
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
    #     Δpᵦₐᵣ  : applied pressure in bar
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

    
    # Pressure correction for applied pressure p_bar /= 0
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

##=======================================================================
##ABETA_HF_DIRI79 <- function (Θᴷ, s, p_bar)
##=======================================================================
const Pᴴᶠᵦ₁ = (
    a₀ =     12.641,
    a₁ = - 1590.2,
    a₂ = -    1.525,
    v₀ = -    9.78,
    v₁ = -    0.0090,
    v₂ = -    0.942E-03,
    k₀ = -    3.91e-3,
    k₁ =      0.054e-3
)
@inline """
    Cᴴᶠᵦ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴᶠᵦ₁)

    # Function calculates association constant \beta_HF [(mol/kg-SW)^-1]
    # in (mol/kg-SW)^-1, where
    #   \beta_HF = \frac [HF]  [H^+] [F^-] 
    # References: Dickson and Riley (1979)
    #             Millero (1995) for pressure correction
    # pH scale  : free
    # Note      : converted here from mol/kg-H2O to mol/kg-SW
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
    #     Δpᵦₐᵣ  : applied pressure in bar
"""
function Fᴴᶠᵦ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴᶠᵦ₁)
    ln_bhf_p0 = ( Pᴴᶠᵦ₁.a₀ +
                  Pᴴᶠᵦ₁.a₁/Θᴷ +
                  Pᴴᶠᵦ₁.a₂*sqrt(μₛ(Sᵖ)))

    # Pressure correction for applied pressure p_bar /= 0
    ln_khf_pp = (-( Pᴴᶠᵦ₁.v₀ +
                     Pᴴᶠᵦ₁.v₁*ΘᵒC(Θᴷ) +
                     Pᴴᶠᵦ₁.v₂*ΘᵒC(Θᴷ)*ΘᵒC(Θᴷ)
                   )+( 
                      Pᴴᶠᵦ₁.k₀ +
                      Pᴴᶠᵦ₁.k₁*ΘᵒC(Θᴷ)
                     )*(Δpᵦₐᵣ/2.)
                     )*(Δpᵦₐᵣ/Rₜ(Θᴷ))
                     
    # Final \beta_HF value
    # --------------------
    #  notice that  ln(k_HF(P)) = ln(k_HF(0)) + zln_khf_pp
    #         <=>  -ln(\beta_HF(P)) = -ln(\beta_HF(0)) + zln_khf_pp
    #         <=>   ln(\beta_HF(P)) =  ln(\beta_HF(0)) - zln_khf_pp

    exp(ln_bhf_p0 - ln_khf_pp ) / H₂Oˢʷ(Sᵖ)
end

##=======================================================================
##AK_HF_PEFR87 <- function (Θᴷ, s, p_bar)
##=======================================================================
const Pᴴᶠₖ₁ = (
    a₀ = -  9.68,
    a₁ =  874.0,
    a₂ =    0.111,
    v₀ = -  9.78,
    v₁ = -  0.0090,
    v₂ = -  0.942E-03,
    k₀ = -  3.91e-3,
    k₁ =    0.054e-3
)
@inline """
    Cᴴᶠₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴᶠₖ₁)

    # Function calculates dissociation constant for hydrogen fluoride
    # in mol/kg-SW
    # References: Perez and Fraga (1987)
    #             Millero (1995) for pressure correction
    # pH scale  : Total (according to Handbook, 2007)
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
    #     Δpᵦₐᵣ  : applied pressure in bar
"""
function Fᴴᶠₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴᶠₖ₁)
    ln_khf_p0 = ( Pᴴᶠₖ₁.a₀ +
                  Pᴴᶠₖ₁.a₁/Θᴷ +
                  Pᴴᶠₖ₁.a₂*Sᴾ⁰⁵(Sᵖ))

    # Pressure correction for applied pressure p_bar /= 0
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

##=======================================================================
##AK_HSO4_DICK90 <- function (Θᴷ, s, p_bar)
##=======================================================================
const Pᴴˢᴼ⁴ₖ₁ = (
    a₀ =     141.328,
    a₁ = -  4276.1,
    a₂ = -    23.093,
    b₀ =     324.57,
    b₁ = - 13856.,
    b₂ = -    47.986,
    c₀ = -   771.54,
    c₁ =   35474.,
    c₂ =     114.723,
    d₀ = -  2698.,
    d₁ =    1776.,
    v₀ = -    18.03,
    v₁ =       0.0466,
    v₂ =       0.316E-03,
    k₀ = -     4.53e-3,
    k₁ =       0.0900e-3
)
@inline """
    Cᴴˢᴼ⁴ₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴˢᴼ⁴ₖ₁)

    # Function returns the dissociation constant of hydrogen sulfate (bisulfate)
    # References: Dickson (1990) -- also Handbook (2007)
    #             Millero (1995) for pressure correction
    # pH scale  : free
    # Note      : converted here from mol/kg-H2O to mol/kg-SW
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
    #     Δpᵦₐᵣ  : applied pressure in bar
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

    # Pressure correction for applied pressure p_bar /= 0
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

##=======================================================================
##ASP_CALC_MUCC83 <- function (Θᴷ, s, p_bar)
##=======================================================================
const Pᶜᵃˡᶜⁱᵗᵉₛₚ = (
    a₀ = - 171.9065,
    a₁ = -   0.077993,
    a₂ =  2839.319,
    a₃ =    71.595,
    b₀ = -   0.77712,
    b₁ =     0.0028426,
    b₂ =   178.34,
    c₀ = -   0.07711,
    d₀ =     0.0041249,
    v₀ = -  48.76,
    v₁ =     0.5304,
    k₀ = -  11.76e-3,
    k₁ =     0.3692e-3
)

@inline """
    Cᶜᵃˡᶜⁱᵗᵉₛₚ(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᶜᵃˡᶜⁱᵗᵉₛₚ)

    # Function returns stoechiometric solubility product
    # of calcite in seawater
    # References: Mucci (1983)
    #             Millero (1995) for pressure correction
    # pH scale  : N/A
    # Units     : (mol/kg-SW)^2
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
    #     Δpᵦₐᵣ  : applied pressure in bar
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

    # Pressure correction for applied pressure p_bar /= 0
    ln_kcalcite_pp = (-( Pᶜᵃˡᶜⁱᵗᵉₛₚ.v₀ +
                            Pᶜᵃˡᶜⁱᵗᵉₛₚ.v₁*ΘᵒC(Θᴷ)
                          )+( 
                             Pᶜᵃˡᶜⁱᵗᵉₛₚ.k₀ +
                             Pᶜᵃˡᶜⁱᵗᵉₛₚ.k₁*ΘᵒC(Θᴷ)
                            )*(Δpᵦₐᵣ/2.)
                            )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

 10^(log10_kcalcite_p0) * exp(ln_kcalcite_pp)
end

##=======================================================================
##ASP_ARAG_MUCC83 <- function (Θᴷ, s, p_bar)
##=======================================================================
const Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ = (
    a₀ = - 171.945,
    a₁ = -   0.077993,
    a₂ =  2903.293,
    a₃ =    71.595,
    b₀ = -   0.068393,
    b₁ =     0.0017276,
    b₂ =    88.135,
    c₀ = -   0.10018,
    d₀ =     0.0059415,
    v₀ = -  48.76,
    v₁ =     0.5304,
    v₂ =     2.8,
    k₀ = -  11.76e-3,
    k₁ =     0.3692e-3
)

@inline """
    Cᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ)

    # Function returns stoechiometric solubility product
    # of aragonite in seawater
    # References: Mucci (1983)
    #             Millero (1979) for pressure correction
    # pH scale  : N/A
    # Units     : (mol/kg-SW)^2
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
    #     Δpᵦₐᵣ  : applied pressure in bar
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

    # Pressure correction for applied pressure p_bar /= 0
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

##=======================================================================
##ACVT_HSWS_O_HTOT <- function (Θᴷ, s, p_bar)
##=======================================================================
@inline """
    H⁺ₛoverH⁺ₜ(Θᴷ,Sᵖ,Δpᵦₐᵣ)

    # Function returns the ratio H_SWS/H_Tot as a function of salinity s
    # Reference:  Munhoven
    # pH scale:   all
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
    #     Δpᵦₐᵣ  : applied pressure in bar
"""
function H⁺ₛoverH⁺ₜ(Θᴷ,Sᵖ,Δpᵦₐᵣ)
    (1. +  
            (Fᵀᴼᵀ(Sᵖ)*Fᴴᶠᵦ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴᶠᵦ₁))
           /(1. + SO₄ᵀᴼᵀ(Sᵖ)/Fᴴˢᴼ⁴ₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴˢᴼ⁴ₖ₁))
           )
end

##=======================================================================
##ACVT_HTOT_O_HFREE <- function (Θᴷ, s, p_bar)
##=======================================================================
@inline """
    H⁺ₜoverH⁺₃(Θᴷ,Sᵖ,Δpᵦₐᵣ)

    # Function returns the ratio H_Tot/H_free as a function of salinity s
    # Reference:  Munhoven
    # pH scale:   N/A
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
    #     Δpᵦₐᵣ  : applied pressure in bar
"""
function H⁺ₜoverH⁺₃(Θᴷ,Sᵖ,Δpᵦₐᵣ)
    1. + SO₄ᵀᴼᵀ(Sᵖ)/Fᴴˢᴼ⁴ₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴˢᴼ⁴ₖ₁)
end

##=======================================================================
##ACVT_HSWS_O_HFREE <- function (Θᴷ, s, p_bar)
##=======================================================================
@inline """
    H⁺ₛoverH⁺₃(Θᴷ,Sᵖ,Δpᵦₐᵣ)

    # Function returns the ratio H_SWS/H_free as a function
    # of salinity s
    # Reference:  Munhoven
    # pH scale:   N/A
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
    #     Δpᵦₐᵣ  : applied pressure in bar
"""
function H⁺ₛoverH⁺₃(Θᴷ,Sᵖ,Δpᵦₐᵣ)
    ( 1. + 
            (Fᵀᴼᵀ(Sᵖ)*Fᴴᶠᵦ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴᶠᵦ₁)) +
            (SO₄ᵀᴼᵀ(Sᵖ)/Fᴴˢᴼ⁴ₖ₁(Θᴷ,Sᵖ,Δpᵦₐᵣ,Pᴴˢᴼ⁴ₖ₁))
         )
end

##=======================================================================
##A_RHOSW1_MUNH97 <- function (Θᴷ, s, p_bar)
##=======================================================================
 const Pᵨ₁ᴹᵘⁿʰ⁹⁷ =
    ( Sᵖ₀ =    35.5,
      Θᴷ₀ =   285.16,
      zₘ₀ =   300.0,
      a₀  =  1039.9044,
      a₁  =     0.77629393,
      a₂  = -   0.19692738
    )

@inline """
    ρ₁ᴹᵘⁿʰ⁹⁷(Θᴷ,Sᵖ,zₘ,Pᵨ₁ᴹᵘⁿʰ⁹⁷)

    # Function returns first order approximation of \rho in (kg-SW)/(m^3-SW)
    # References: Munhoven (1997)
    #             after EOS80 (UNESCO, 1981, 1983)
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
"""
function ρ₁ᴹᵘⁿʰ⁹⁷(Θᴷ,Sᵖ,zₘ,Pᵨ₁ᴹᵘⁿʰ⁹⁷)
    ( Pᵨ₁ᴹᵘⁿʰ⁹⁷.a₀ + 
             Pᵨ₁ᴹᵘⁿʰ⁹⁷.a₁*(Sᵖ-Pᵨ₁ᴹᵘⁿʰ⁹⁷.Sᵖ₀) +
             Pᵨ₁ᴹᵘⁿʰ⁹⁷.a₂*(Θᴷ-Pᵨ₁ᴹᵘⁿʰ⁹⁷.Θᴷ0) 
           )
end

##=======================================================================
##A_RHOSW2_MUNH97 <- function (Θᴷ, s, p_bar)
##=======================================================================
const Pᵨ₂ᴹᵘⁿʰ⁹⁷ =
    ( Sᵖ₀ =     35.5,
      Θᴷ₀ =    285.16,
      zₘ₀ =    300.0,
      ρ₀  =   1040.0145,
      ρ₁  =      0.77629393,
      ρ₂  = -    0.25013591,
      ρ₃  =      0.042026266e-2,
      ρ₄  = -    0.047473116e-3,
      ρ₅  = -    0.047974224e-6,
      ρ₆  = -    0.21404592e-4
    )

@inline """
    ρᴹᵘⁿʰ⁹⁷(Θᴷ,Sᵖ,zₘ,Pᵨ₂ᴹᵘⁿʰ⁹⁷)

    # Function returns second order approximation of \rho in (kg-SW)/(m^3-SW)
    # References: Munhoven (1997)
    #             after EOS80 (UNESCO, 1981, 1983)
    # ------------------
    # Argument variables
    # ------------------
    #     Sᵖ    : practical salinity
    #     Θᴷ    : temperature in K
    #     zₘ    : depth in metres
"""
function ρᴹᵘⁿʰ⁹⁷(Θᴷ,Sᵖ,zₘ,Pᵨ₂ᴹᵘⁿʰ⁹⁷)
    (
        Pᵨ₂ᴹᵘⁿʰ⁹⁷.ρ₀ + 
        Pᵨ₂ᴹᵘⁿʰ⁹⁷.ρ₁*(Sᵖ-Pᵨ₂ᴹᵘⁿʰ⁹⁷.Sᵖ₀) +
        Pᵨ₂ᴹᵘⁿʰ⁹⁷.ρ₂*(Θᴷ-Pᵨ₂ᴹᵘⁿʰ⁹⁷.Θᴷ₀) +
        Pᵨ₂ᴹᵘⁿʰ⁹⁷.ρ₃*(zₘ-Pᵨ₂ᴹᵘⁿʰ⁹⁷.zₘ₀) +
        Pᵨ₂ᴹᵘⁿʰ⁹⁷.ρ₄*(Θᴷ-Pᵨ₂ᴹᵘⁿʰ⁹⁷.Θᴷ₀)*(Θᴷ-Pᵨ₂ᴹᵘⁿʰ⁹⁷.Θᴷ₀) +
        Pᵨ₂ᴹᵘⁿʰ⁹⁷.ρ₅*(zₘ-Pᵨ₂ᴹᵘⁿʰ⁹⁷.zₘ₀)*(zₘ-Pᵨ₂ᴹᵘⁿʰ⁹⁷.zₘ₀) +
        Pᵨ₂ᴹᵘⁿʰ⁹⁷.ρ₆*(Θᴷ-Pᵨ₂ᴹᵘⁿʰ⁹⁷.Θᴷ₀)*(zₘ-Pᵨ₂ᴹᵘⁿʰ⁹⁷.zₘ₀)
    )
end

end # module DissociationConstants
#------------------------------------------------------------
module CarbonSystemApprox
export CarbonSolverApprox

using ..DissociationConstants
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
using .DissociationConstants

## This should go in the testing suite, eventually.
Θᶜ      = 25.0
Sᴬ      = 35.0
Δpᵦₐᵣ   = 0.0
Cᵀ      = 2050.0*1e-6 # umol/kg to mol/kg
Aᵀ      = 2350.0*1e-6 # umol/kg to mol/kg
pCO₂ᵃᵗᵐ = 280.0*1e-6  # uatm to atm
pH      = 8.0

Cᶜᵒⁿˢᵗ = DissociationCoefficients(Θᶜ, Sᴬ, Δpᵦₐᵣ);

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

pCO₂, pH = CarbonSolverApprox(
        Θᶜ, Sᴬ, Δpᵦₐᵣ, Cᵀ, Aᵀ, pH, pCO₂ᵃᵗᵐ,
        )
println("Cᵀ = ", Cᵀ*1e6    )
println("Aᵀ = ", Aᵀ*1e6    )
println("pH = "  , pH      )
println("pCO₂ᵃᵗᵐ = ", pCO₂ᵃᵗᵐ*1e6)
println("pCO₂ = ", pCO₂*1e6)

end # module CarbonSolver