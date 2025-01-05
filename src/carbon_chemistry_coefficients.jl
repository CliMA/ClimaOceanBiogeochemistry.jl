struct CarbonChemistryCoefficients{FT}
    Cᵈⁱᶜₖₛₒₗₐ     :: FT
    Cᵈⁱᶜₖₛₒₗₒ     :: FT
    Cᵈⁱᶜₖ₀       :: FT 
    Cᵈⁱᶜₖ₁ᵣ₉₃    :: FT 
    Cᵈⁱᶜₖ₂ᵣ₉₃    :: FT 
    Cᵈⁱᶜₖ₁ₘ₉₅    :: FT 
    Cᵈⁱᶜₖ₂ₘ₉₅    :: FT 
    Cᵈⁱᶜₖ₁ₗ₀₀     :: FT
    Cᵈⁱᶜₖ₂ₗ₀₀     :: FT
    Cᵇₖ₁         :: FT 
    Cᴴ²ᴼₖ₁       :: FT 
    Cᴾᴼ⁴ₖ₁       :: FT 
    Cᴾᴼ⁴ₖ₂       :: FT 
    Cᴾᴼ⁴ₖ₃       :: FT 
    Cˢⁱᵗₖ₁       :: FT 
    Cᴴ²ˢₖ₁       :: FT 
    Cᴺᴴ⁴ₖ₁       :: FT
    Cᴴᶠᵦ₁        :: FT
    Cᴴᶠₖ₁        :: FT
    Cᴴˢᴼ⁴ₖ₁      :: FT
    Cᶜᵃˡᶜⁱᵗᵉₛₚ   :: FT
    Cᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ :: FT
    Cᴮᵀ         :: FT 
    Cᶠᵀ         :: FT 
    Cᶜᵃ         :: FT 
    Cˢᴼ⁴        :: FT 
    H⁺ₛoverH⁺ₜ   :: FT
    H⁺ₜoverH⁺₃   :: FT
    H⁺ₛoverH⁺₃   :: FT
end

"""
    CarbonChemistryCoefficients(Θᶜ, Sᴬ, Δpᵦₐᵣ, params)

Return dissociation coefficients necessary to solve for the distribution of carbonate species.
"""
@inline function CarbonChemistryCoefficients(
    params :: CarbonSystemParameters,
    Θᶜ     :: Real = 25, 
    Sᴬ     :: Real = 35, 
    Δpᵦₐᵣ  :: Real = 0, 
    )
    
# Need a conversion from Absolute Salinity, Sᴬ, to Practical Salinity, Sᴾ
    Sᵖ = Sᴬ
# Sᴾ = Sᴬ / (1.0 - 0.35 * Sᴬ / 35.0) ??
# What about converting temperature from Conservative Temperature to potential temperature?
    Θᴷ = Θᵒᴷ(Θᶜ)
# Also need to convert from Absolute Pressure, Pᴬ, to Applied Pressure in bars, the pressure relative to (1x) atm

# Calculate the coefficients
    Cᵈⁱᶜₖₛₒₗₐ    = Fᵈⁱᶜₖₛₒₗₐ(Θᴷ, Sᵖ, params.Pᵈⁱᶜₖₛₒₗₐ)
    Cᵈⁱᶜₖₛₒₗₒ    = Fᵈⁱᶜₖₚᵣₑ(Θᴷ, Sᵖ, params.Pᵈⁱᶜₖₚᵣₑ) * Fᵈⁱᶜₖ₀(Θᴷ, Sᵖ, params.Pᵈⁱᶜₖ₀) # Fᵈⁱᶜₖₛₒₗₒ
    Cᵈⁱᶜₖ₀       = Fᵈⁱᶜₖ₀(Θᴷ, Sᵖ, params.Pᵈⁱᶜₖ₀)
    Cᵈⁱᶜₖ₁ᵣ₉₃    = Fᵈⁱᶜₖ₁ᵣ₉₃(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᵈⁱᶜₖ₁ᵣ₉₃, params.Pᴴ²⁰ˢʷ)
    Cᵈⁱᶜₖ₂ᵣ₉₃    = Fᵈⁱᶜₖ₂ᵣ₉₃(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᵈⁱᶜₖ₂ᵣ₉₃, params.Pᴴ²⁰ˢʷ)
    Cᵈⁱᶜₖ₁ₘ₉₅    = Fᵈⁱᶜₖ₁ₘ₉₅(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᵈⁱᶜₖ₁ₘ₉₅)/H⁺ₛoverH⁺ₜ(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᶠᵀᴼᵀ, params.Pᴴᶠᵦ₁, params.Pᵘˢ, params.Pᴴ²⁰ˢʷ, params.Pˢᴼ⁴ᵀᴼᵀ, params.Pᴴˢᴼ⁴ₖ₁)
    Cᵈⁱᶜₖ₂ₘ₉₅    = Fᵈⁱᶜₖ₂ₘ₉₅(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᵈⁱᶜₖ₂ₘ₉₅)/H⁺ₛoverH⁺ₜ(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᶠᵀᴼᵀ, params.Pᴴᶠᵦ₁, params.Pᵘˢ, params.Pᴴ²⁰ˢʷ, params.Pˢᴼ⁴ᵀᴼᵀ, params.Pᴴˢᴼ⁴ₖ₁)
    Cᵈⁱᶜₖ₁ₗ₀₀    = Fᵈⁱᶜₖ₁ₗ₀₀(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᵈⁱᶜₖ₁ₗ₀₀)
    Cᵈⁱᶜₖ₂ₗ₀₀    = Fᵈⁱᶜₖ₂ₗ₀₀(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᵈⁱᶜₖ₂ₗ₀₀)
    Cᵇₖ₁         = Fᵇₖ₁(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᴮₖ₁)
    Cᴴ²ᴼₖ₁       = Fᴴ²ᴼₖ₁(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᴴ²ᴼₖ₁)/H⁺ₛoverH⁺ₜ(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᶠᵀᴼᵀ, params.Pᴴᶠᵦ₁, params.Pᵘˢ, params.Pᴴ²⁰ˢʷ, params.Pˢᴼ⁴ᵀᴼᵀ, params.Pᴴˢᴼ⁴ₖ₁)
    Cᴾᴼ⁴ₖ₁       = Fᴾᴼ⁴ₖ₁(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᴾᴼ⁴ₖ₁)/H⁺ₛoverH⁺ₜ(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᶠᵀᴼᵀ, params.Pᴴᶠᵦ₁, params.Pᵘˢ, params.Pᴴ²⁰ˢʷ, params.Pˢᴼ⁴ᵀᴼᵀ, params.Pᴴˢᴼ⁴ₖ₁)
    Cᴾᴼ⁴ₖ₂       = Fᴾᴼ⁴ₖ₂(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᴾᴼ⁴ₖ₂)/H⁺ₛoverH⁺ₜ(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᶠᵀᴼᵀ, params.Pᴴᶠᵦ₁, params.Pᵘˢ, params.Pᴴ²⁰ˢʷ, params.Pˢᴼ⁴ᵀᴼᵀ, params.Pᴴˢᴼ⁴ₖ₁)
    Cᴾᴼ⁴ₖ₃       = Fᴾᴼ⁴ₖ₃(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᴾᴼ⁴ₖ₃)/H⁺ₛoverH⁺ₜ(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᶠᵀᴼᵀ, params.Pᴴᶠᵦ₁, params.Pᵘˢ, params.Pᴴ²⁰ˢʷ, params.Pˢᴼ⁴ᵀᴼᵀ, params.Pᴴˢᴼ⁴ₖ₁)
    Cˢⁱᵗₖ₁       = Fˢⁱᵗₖ₁(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pˢⁱᵗₖ₁, params.Pᵘˢ, params.Pᴴ²⁰ˢʷ)/H⁺ₛoverH⁺ₜ(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᶠᵀᴼᵀ, params.Pᴴᶠᵦ₁, params.Pᵘˢ, params.Pᴴ²⁰ˢʷ, params.Pˢᴼ⁴ᵀᴼᵀ, params.Pᴴˢᴼ⁴ₖ₁)
    Cᴴ²ˢₖ₁       = Fᴴ²ˢₖ₁(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᴴ²ˢₖ₁)/H⁺ₛoverH⁺ₜ(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᶠᵀᴼᵀ, params.Pᴴᶠᵦ₁, params.Pᵘˢ, params.Pᴴ²⁰ˢʷ, params.Pˢᴼ⁴ᵀᴼᵀ, params.Pᴴˢᴼ⁴ₖ₁)
    Cᴺᴴ⁴ₖ₁       = Fᴺᴴ⁴ₖ₁(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᴺᴴ⁴ₖ₁)/H⁺ₛoverH⁺ₜ(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᶠᵀᴼᵀ, params.Pᴴᶠᵦ₁, params.Pᵘˢ, params.Pᴴ²⁰ˢʷ, params.Pˢᴼ⁴ᵀᴼᵀ, params.Pᴴˢᴼ⁴ₖ₁)
    Cᴴᶠᵦ₁        = Fᴴᶠᵦ₁(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᴴᶠᵦ₁, params.Pᵘˢ, params.Pᴴ²⁰ˢʷ)
    Cᴴᶠₖ₁        = Fᴴᶠₖ₁(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᴴᶠₖ₁)
    Cᴴˢᴼ⁴ₖ₁      = Fᴴˢᴼ⁴ₖ₁(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᴴˢᴼ⁴ₖ₁, params.Pᵘˢ, params.Pᴴ²⁰ˢʷ) # Leave on the free scale
    Cᶜᵃˡᶜⁱᵗᵉₛₚ   = Fᶜᵃˡᶜⁱᵗᵉₛₚ(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᶜᵃˡᶜⁱᵗᵉₛₚ)
    Cᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ = Fᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ)
    Cᴮᵀ          = Bᵀᴼᵀ(Sᵖ, params.Pᴮᵀᴼᵀ)
    Cᶠᵀ          = Fᵀᴼᵀ(Sᵖ, params.Pᶠᵀᴼᵀ)
    Cᶜᵃ          = Caᵀᴼᵀ(Sᵖ, params.Pᶜᵃᵀᴼᵀ)
    Cˢᴼ⁴         = SO₄ᵀᴼᵀ(Sᵖ, params.Pˢᴼ⁴ᵀᴼᵀ)
    CH⁺ₛoverH⁺ₜ   = H⁺ₛoverH⁺ₜ(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᶠᵀᴼᵀ, params.Pᴴᶠᵦ₁, params.Pᵘˢ, params.Pᴴ²⁰ˢʷ, params.Pˢᴼ⁴ᵀᴼᵀ, params.Pᴴˢᴼ⁴ₖ₁)
    CH⁺ₜoverH⁺₃   = H⁺ₜoverH⁺₃(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pˢᴼ⁴ᵀᴼᵀ, params.Pᴴˢᴼ⁴ₖ₁, params.Pᵘˢ, params.Pᴴ²⁰ˢʷ)
    CH⁺ₛoverH⁺₃   = H⁺ₛoverH⁺₃(Θᴷ, Sᵖ, Δpᵦₐᵣ, params.Pᶠᵀᴼᵀ, params.Pᴴᶠᵦ₁, params.Pᵘˢ, params.Pᴴ²⁰ˢʷ, params.Pˢᴼ⁴ᵀᴼᵀ, params.Pᴴˢᴼ⁴ₖ₁)

    return CarbonChemistryCoefficients(
        Cᵈⁱᶜₖₛₒₗₐ,
        Cᵈⁱᶜₖₛₒₗₒ,  
        Cᵈⁱᶜₖ₀,  
        Cᵈⁱᶜₖ₁ᵣ₉₃,   
        Cᵈⁱᶜₖ₂ᵣ₉₃,   
        Cᵈⁱᶜₖ₁ₘ₉₅,   
        Cᵈⁱᶜₖ₂ₘ₉₅,   
        Cᵈⁱᶜₖ₁ₗ₀₀,   
        Cᵈⁱᶜₖ₂ₗ₀₀,   
        Cᵇₖ₁,        
        Cᴴ²ᴼₖ₁,      
        Cᴾᴼ⁴ₖ₁,      
        Cᴾᴼ⁴ₖ₂,      
        Cᴾᴼ⁴ₖ₃,      
        Cˢⁱᵗₖ₁,      
        Cᴴ²ˢₖ₁,      
        Cᴺᴴ⁴ₖ₁,      
        Cᴴᶠᵦ₁,      
        Cᴴᶠₖ₁,       
        Cᴴˢᴼ⁴ₖ₁,     
        Cᶜᵃˡᶜⁱᵗᵉₛₚ,  
        Cᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ,
        Cᴮᵀ,         
        Cᶠᵀ,         
        Cᶜᵃ,         
        Cˢᴼ⁴,        
        CH⁺ₛoverH⁺ₜ,  
        CH⁺ₜoverH⁺₃,  
        CH⁺ₛoverH⁺₃,  
       )
end

adapt_structure( 
        to, ccc::CarbonChemistryCoefficients,
    ) = CarbonChemistryCoefficients(
    adapt(to, ccc.Cᵈⁱᶜₖₛₒₗₐ),
    adapt(to, ccc.Cᵈⁱᶜₖₛₒₗₒ),
    adapt(to, ccc.Cᵈⁱᶜₖ₀),
    adapt(to, ccc.Cᵈⁱᶜₖ₁ᵣ₉₃),
    adapt(to, ccc.Cᵈⁱᶜₖ₂ᵣ₉₃),
    adapt(to, ccc.Cᵈⁱᶜₖ₁ₘ₉₅),
    adapt(to, ccc.Cᵈⁱᶜₖ₂ₘ₉₅),
    adapt(to, ccc.Cᵈⁱᶜₖ₁ₗ₀₀),
    adapt(to, ccc.Cᵈⁱᶜₖ₂ₗ₀₀),
    adapt(to, ccc.Cᵇₖ₁),
    adapt(to, ccc.Cᴴ²ᴼₖ₁),
    adapt(to, ccc.Cᴾᴼ⁴ₖ₁),
    adapt(to, ccc.Cᴾᴼ⁴ₖ₂),
    adapt(to, ccc.Cᴾᴼ⁴ₖ₃),
    adapt(to, ccc.Cˢⁱᵗₖ₁),
    adapt(to, ccc.Cᴴ²ˢₖ₁),
    adapt(to, ccc.Cᴺᴴ⁴ₖ₁),
    adapt(to, ccc.Cᴴᶠᵦ₁),
    adapt(to, ccc.Cᴴᶠₖ₁),
    adapt(to, ccc.Cᴴˢᴼ⁴ₖ₁),
    adapt(to, ccc.Cᶜᵃˡᶜⁱᵗᵉₛₚ),
    adapt(to, ccc.Cᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ),
    adapt(to, ccc.Cᴮᵀ),
    adapt(to, ccc.Cᶠᵀ),
    adapt(to, ccc.Cᶜᵃ),
    adapt(to, ccc.Cˢᴼ⁴),
    adapt(to, ccc.H⁺ₛoverH⁺ₜ),
    adapt(to, ccc.H⁺ₜoverH⁺₃),
    adapt(to, ccc.H⁺ₛoverH⁺₃),
)

# Gas constant
const gasconst_bar_cm3_o_mol_k = 83.14472 # Handbook (2007)

# 0 degrees centigrade in Kelvin
const Θᴷ_0ᵒC = 273.15 # Handbook (2007)

# Salinity reference value
const Sʳᵉᶠ = -34.8

@inline Θᵒᴷ(ΘᵒC::Real)::Real = ΘᵒC + Θᴷ_0ᵒC
@inline Θᴷ₁₀₀(Θᴷ::Real)::Real = Θᴷ/100
@inline Sᴾ⁰⁵(Sᵖ::Real)::Real  = sqrt(Sᵖ)
@inline ΘᵒC(Θᴷ::Real)::Real   = Θᴷ - Θᴷ_0ᵒC
@inline ΔSᵖ(Sᵖ::Real)::Real = Sᵖ + Sʳᵉᶠ
@inline Rₜ(Θᴷ::Real)::Real = gasconst_bar_cm3_o_mol_k * Θᴷ

#Base.@kwdef struct Pᴴ²⁰ˢʷ{FT}
#    a₀ :: FT    =   1.0
#    a₁ :: FT    = - 0.001005
#end
#adapt_structure( 
#        to, p::Pᴴ²⁰ˢʷ,
#    ) = Pᴴ²⁰ˢʷ(
#    adapt(to, p.a₀),
#    adapt(to, p.a₁),
#)
"""
    H₂Oˢʷ(Sᵖ, Pᴴ²⁰ˢʷ)

Return the mass of pure water in one kg of seawater
of practical salinity, `Sᵖ`.
References: "libthdyct" -- derived by Munhoven (1997) from data by Millero (1982)
            "Handbook (2007)" -- Handbook (2007)
pH scale:   N/A
"""
@inline function H₂Oˢʷ(Sᵖ::Real, p₁ = Pᴴ²⁰ˢʷ) :: Real
    (; a₀, a₁) = p₁
    return a₀ + a₁ * Sᵖ 
end
#H₂Oˢʷ(Sᵖ, Pᴴ²⁰ˢʷ) = 1. - 0.0010049*Sᵖ # libthdyct

# Base.@kwdef struct Pᵘˢ{FT}
#     a₀ :: FT    =   0.019924
# end
# adapt_structure( 
#         to, p::Pᵘˢ,
#     ) = Pᵘˢ(
#     adapt(to, p.a₀),
# )
"""
    μₛ(Sᵖ)

Return ionic strength in mol/kg-SW, for given practical salinity, `Sᵖ`.
References: "libthdyct" -- derived by Munhoven (1997) from data by Millero (1982)
            "Handbook (2007)" -- Handbook (2007)
pH scale:   N/A
"""
@inline μₛ(Sᵖ::Real, p₁ = Pᵘˢ, p₂ = Pᴴ²⁰ˢʷ) ::Real = (p₁.a₀ * Sᵖ) / H₂Oˢʷ(Sᵖ, p₂) # Handbook (2007)
# μₛ(Sᵖ)    = (0.019920 * Sᵖ) / H₂Oˢʷ(Sᵖ, Pᴴ²⁰ˢʷ)# libthdyct

# Base.@kwdef struct Pᴮᵀᴼᵀ{FT}
#     a₀ :: FT    =   0.000416
#     a₁ :: FT    =   35.0
#     a₂ :: FT    =   1.0
# end
# adapt_structure( 
#         to, p::Pᴮᵀᴼᵀ,
#     ) = Pᴮᵀᴼᵀ(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
# )
"""
    Bᵀᴼᵀ(Sᵖ, Pᴮᵀᴼᵀ)

Return total borate concentration in mol/kg-SW given practical salinity, `Sᵖ`.
References: Uppström (1974), cited by  Dickson et al. (2007, chapter 5, p 10)
            Millero (1982) cited in Millero (1995)
"""
@inline function Bᵀᴼᵀ(Sᵖ::Real, p₁ = Pᴮᵀᴼᵀ) :: Real
    (; a₀, a₁, a₂) = p₁
    return a₀ * (Sᵖ / a₁) / a₂
end
#Bᵀᴼᵀ(Sᵖ, Pᴮᵀᴼᵀ)  = 0.000232 * (Sᵖ/1.80655)/10.811

# Base.@kwdef struct Pᶜᵃᵀᴼᵀ{FT}
#     a₀ :: FT    =   0.02127
#     a₁ :: FT    =   40.078
#     a₂ :: FT    =   1.80655
# end
# adapt_structure( 
#         to, p::Pᶜᵃᵀᴼᵀ,
#     ) = Pᶜᵃᵀᴼᵀ(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
# )
"""
    Caᵀᴼᵀ(Sᵖ, Pᶜᵃᵀᴼᵀ)

Return calcium concentration in mol/kg-SW given practical salinity, `Sᵖ`.
References: Culkin and Cox (1966), 
            Culkin (1967), 
            Riley and Tongudai (1967)    
"""
@inline function Caᵀᴼᵀ(Sᵖ::Real, p₁ = Pᶜᵃᵀᴼᵀ) :: Real
    (; a₀, a₁, a₂) = p₁
    return (a₀ / a₁) * (Sᵖ / a₂)
end
# Caᵀᴼᵀ(Sᵖ, Pᶜᵃᵀᴼᵀ) = 0.010282*(Sᵖ/35.)

# Base.@kwdef struct Pᶠᵀᴼᵀ{FT}
#     a₀ :: FT    =   6.8e-5
#     a₁ :: FT    =   35.0
# end
# adapt_structure( 
#         to, p::Pᶠᵀᴼᵀ,
#     ) = Pᶠᵀᴼᵀ(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
# )
"""
    Fᵀᴼᵀ(Sᵖ, Pᶠᵀᴼᵀ)

Return total fluoride concentration in mol/kg-SW given practical salinity, `Sᵖ`.
References: Culkin (1965) (???)
"""
@inline function Fᵀᴼᵀ(Sᵖ::Real, p₁ = Pᶠᵀᴼᵀ) :: Real 
    (; a₀, a₁) = p₁
    return a₀ * (Sᵖ / a₁)
end

# Base.@kwdef struct Pˢᴼ⁴ᵀᴼᵀ{FT}
#     a₀ :: FT    =   0.1400
#     a₁ :: FT    =   96.062
#     a₂ :: FT    =   1.80655
# end
# adapt_structure( 
#         to, p::Pˢᴼ⁴ᵀᴼᵀ,
#     ) = Pˢᴼ⁴ᵀᴼᵀ(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
# )
"""
    SO₄ᵀᴼᵀ(Sᵖ, Pˢᴼ⁴ᵀᴼᵀ)

Return total sulfate concentration in mol/kg-SW given practical salinity, `Sᵖ`.
References: Morris, A.W. and Riley, J.P. (1966) quoted in Handbook (2007)
"""
@inline function SO₄ᵀᴼᵀ(Sᵖ::Real, p₁ = Pˢᴼ⁴ᵀᴼᵀ) :: Real
    (; a₀, a₁, a₂) = p₁
    return (a₀ / a₁) * (Sᵖ / a₂)
end
#SO₄ᵀᴼᵀ(Sᵖ, Pˢᴼ⁴ᵀᴼᵀ) = 0.028234*(Sᵖ/35.)

#Base.@kwdef struct Pᴺᴴ⁴ᵀᴼᵀ{FT}
#    a₀ :: FT    =   0.0001
#    a₁ :: FT    =   17.031
#    a₂ :: FT    =   1.80655
#end
# I cannot independently verify this equation.  It is not in the Handbook (2007)
#"""
#NH₄ᵀᴼᵀ(Sᵖ, Pᴺᴴ⁴ᵀᴼᵀ)
#
#Return total ammonium concentration in mol/kg-SW given practical salinity, `Sᵖ`.
#References: Yamamoto (1995)
#"""
#@inline function NH₄ᵀᴼᵀ(Sᵖ, p₁ = Pᴺᴴ⁴ᵀᴼᵀ) 
#    (; a₀, a₁, a₂) = p₁()
#    return (a₀ / a₁) * (Sᵖ / a₂)
#end
#
#Base.@kwdef struct Pᴾᴼ⁴ᵀᴼᵀ{FT}
#    a₀ :: FT    =   0.0001
#    a₁ :: FT    =   94.971
#    a₂ :: FT    =   1.80655
#end
#
#"""
#PO₄ᵀᴼᵀ(Sᵖ, Pᴾᴼ⁴ᵀᴼᵀ)
#
#Return total phosphate concentration in mol/kg-SW given practical salinity, `Sᵖ`.
#References: Millero (1995)
#"""
#@inline function PO₄ᵀᴼᵀ(Sᵖ, p₁ = Pᴾᴼ⁴ᵀᴼᵀ) 
#    (; a₀, a₁, a₂) = p₁()
#    return (a₀ / a₁) * (Sᵖ / a₂)
#end
#
#Base.@kwdef struct Pˢⁱᴼ³ᵀᴼᵀ{FT}
#    a₀ :: FT    =   0.0001
#    a₁ :: FT    =   60.084
#    a₂ :: FT    =   1.80655
#end
#
#"""
#SiO₃ᵀᴼᵀ(Sᵖ, Pˢⁱᴼ³ᵀᴼᵀ)
#
#Return total silicate concentration in mol/kg-SW given practical salinity, `Sᵖ`.
#References: Millero (1995)
#"""
#@inline function SiO₃ᵀᴼᵀ(Sᵖ, p₁ = Pˢⁱᴼ³ᵀᴼᵀ) 
#    (; a₀, a₁, a₂) = p₁()
#    return (a₀ / a₁) * (Sᵖ / a₂)
#end
#
#Base.@kwdef struct Pᴴ²ˢᵀᴼᵀ{FT}
#    a₀ :: FT    =   0.0001
#    a₁ :: FT    =   34.082
#    a₂ :: FT    =   1.80655
#end
#
#"""
#H₂Sᵀᴼᵀ(Sᵖ, Pᴴ²ˢᵀᴼᵀ)
#
#Return total hydrogen sulfide concentration in mol/kg-SW given practical salinity, `Sᵖ`.
#References: Dickson (1990)
#"""
#@inline function H₂Sᵀᴼᵀ(Sᵖ, p₁ = Pᴴ²ˢᵀᴼᵀ) 
#    (; a₀, a₁, a₂) = p₁()
#    return (a₀ / a₁) * (Sᵖ / a₂)
#end

# Base.@kwdef struct Pᵈⁱᶜₖₛₒₗₐ{FT}
#     a₀ :: FT = -162.8301
#     a₁ :: FT =  218.2968
#     a₂ :: FT =   90.9241
#     a₃ :: FT = -  1.47696
#     b₀ :: FT =   0.025695
#     b₁ :: FT = - 0.025225
#     b₂ :: FT =   0.0049867
# end
# adapt_structure( 
#         to, p::Pᵈⁱᶜₖₛₒₗₐ,
#     ) = Pᵈⁱᶜₖₛₒₗₐ(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
#     adapt(to, p.a₃),
#     adapt(to, p.b₀),
#     adapt(to, p.b₁),
#     adapt(to, p.b₂),
# )
"""
    Fᵈⁱᶜₖₛₒₗₐ(Θᴷ, Sᵖ, Pᵈⁱᶜₖₛₒₗₐ)

    Calculate f = k0(1-pH2O)*correction term for non-ideality in (mol/kg-SW)/atm given temperature 
    in K, `Θᴷ`, practical salinity, `Sᵖ`, and coefficients, `Pᵈⁱᶜₖₛₒₗₐ`. Currently no pressure correction

    References: Weiss & Price (1980, Mar. Chem., 8, 347-359 Eq 13 with table 6 values)
    pH scale  : N/A
    Note      : currently no pressure correction
"""
@inline function Fᵈⁱᶜₖₛₒₗₐ(Θᴷ::Real, Sᵖ::Real, p₁ = Pᵈⁱᶜₖₛₒₗₐ) :: Real
    (; a₀, a₁, a₂, a₃, b₀, b₁, b₂) = p₁
    return exp(
               a₀ + 
               a₁/Θᴷ₁₀₀(Θᴷ) +
               a₂*log(Θᴷ₁₀₀(Θᴷ)) +
               a₃*Θᴷ₁₀₀(Θᴷ)*Θᴷ₁₀₀(Θᴷ) +    
               (
                b₀ + 
                b₁*(Θᴷ₁₀₀(Θᴷ)) +
                b₂*Θᴷ₁₀₀(Θᴷ)*Θᴷ₁₀₀(Θᴷ)
              )*Sᵖ
             )
end

# Base.@kwdef struct Pᵈⁱᶜₖₚᵣₑ{FT}
#     a₀ :: FT = -1636.75
#     a₁ :: FT = -  12.0408
#     a₂ :: FT = -   0.0327957 
#     a₃ :: FT =     3.16528e-5
#     b₀ :: FT =    57.7
#     b₁ :: FT = -   0.118
#     p₀ :: FT =     1.01325 # p_bar_oneatmosphere, Handbook (2007)
# end
# adapt_structure( 
#         to, p::Pᵈⁱᶜₖₚᵣₑ,
#     ) = Pᵈⁱᶜₖₚᵣₑ(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
#     adapt(to, p.a₃),
#     adapt(to, p.b₀),
#     adapt(to, p.b₁),
#     adapt(to, p.p₀),
# )
"""
    Fᵈⁱᶜₖₚᵣₑ(Θᴷ, Sᵖ, Pᵈⁱᶜₖₚᵣₑ)

Return fugacity prefactor needed for non-ideality of CO₂ in the  ocean 
    in (mol/kg-SW)/atm given temperature in K, `Θᴷ`, practical salinity, 
    `Sᵖ`, and coefficients, `Pᵈⁱᶜₖₚᵣₑ`.

References: Weiss (1974) Marine Chemistry
pH scale  : N/A
Note      : currently no pressure correction
"""
@inline function Fᵈⁱᶜₖₚᵣₑ(Θᴷ::Real, Sᵖ, p₁ = Pᵈⁱᶜₖₚᵣₑ) :: Real
    (; a₀, a₁, a₂, a₃, b₀, b₁, p₀) = p₁

#  "x2" term often neglected (assumed=1) in applications of Weiss (1974) eq.9
#   x2 = 1 - x1 = 1 - xCO2 (it is very close to 1, but not quite)
    return exp( 
        (a₀ + 
         a₁*Θᴷ + 
         a₂*Θᴷ*Θᴷ + 
         a₃*Θᴷ*Θᴷ*Θᴷ+
        2*(b₀ + b₁*Θᴷ)) * 
        (p₀ / Rₜ(Θᴷ)))
end

# Base.@kwdef struct Pᵈⁱᶜₖ₀{FT}
#     a₀ :: FT = -60.2409
#     a₁ :: FT =  93.4517
#     a₂ :: FT =  23.3585
#     b₀ :: FT =   0.023517
#     b₁ :: FT = - 0.023656
#     b₂ :: FT =   0.0047036
# end
# adapt_structure( 
#         to, p::Pᵈⁱᶜₖ₀,
#     ) = Pᵈⁱᶜₖ₀(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
#     adapt(to, p.b₀),
#     adapt(to, p.b₁),
#     adapt(to, p.b₂),
# )
"""
    Fᵈⁱᶜₖ₀(Θᴷ, Sᵖ, Pᵈⁱᶜₖ₀)

Return hydration constant of CO₂ in (mol/kg-SW)/atm given temperature 
in K, `Θᴷ`, practical salinity, `Sᵖ`, and coefficients, `Pᵈⁱᶜₖ₀`.

CO₂ + H₂O <-> H₂CO₃

References: Weiss (1979)
pH scale  : N/A
Note      : currently no pressure correction
"""
@inline function Fᵈⁱᶜₖ₀(Θᴷ::Real, Sᵖ::Real, p₁ = Pᵈⁱᶜₖ₀) :: Real
    (; a₀, a₁, a₂, b₀, b₁, b₂) = p₁
    return exp(
               a₀ + 
               a₁/Θᴷ₁₀₀(Θᴷ) +
               a₂*log(Θᴷ₁₀₀(Θᴷ)) +    
               (
                b₀ + 
                b₁*(Θᴷ₁₀₀(Θᴷ)) +
                b₂*Θᴷ₁₀₀(Θᴷ)*Θᴷ₁₀₀(Θᴷ)
              )*Sᵖ
             )
end

# Base.@kwdef struct Pᵈⁱᶜₖ₁ᵣ₉₃{FT}
#     a₀ :: FT =     2.83655
#     a₁ :: FT = -2307.1266
#     a₂ :: FT = -   1.5529413
#     b₀ :: FT = -   0.20760841
#     b₁ :: FT = -   4.0484
#     b₂ :: FT =     0.08468345
#     b₃ :: FT = -   0.00654208
#     v₀ :: FT = -  25.5
#     v₁ :: FT = -   0.151
#     v₂ :: FT =     0.1271
#     k₀ :: FT = -   3.08e-3
#     k₁ :: FT = -   0.578e-3
#     k₂ :: FT =     0.0877e-3
# end
# adapt_structure( 
#         to, p::Pᵈⁱᶜₖ₁ᵣ₉₃,
#     ) = Pᵈⁱᶜₖ₁ᵣ₉₃(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
#     adapt(to, p.b₀),
#     adapt(to, p.b₁),
#     adapt(to, p.b₂),
#     adapt(to, p.b₃),
#     adapt(to, p.v₀),
#     adapt(to, p.v₁),
#     adapt(to, p.v₂),
#     adapt(to, p.k₀),
#     adapt(to, p.k₁),
#     adapt(to, p.k₂),
# )
"""
    Fᵈⁱᶜₖ₁ᵣ₉₃(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, Pᵈⁱᶜₖ₁ᵣ₉₃)

Return the first dissociation constant of carbonic acid in mol/kg-SW, given
temperature in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, 
`Δpᵦₐᵣ`, and coefficients, `Pᵈⁱᶜₖ₁ᵣ₉₃`.

H₂CO₃ <-> HCO₃⁻ + H⁺

References: Roy et al. (1993) -- also Handbook (1994)
            Millero (1979) pressure correction
pH scale  : Total
Valid range: T:  0-45  S:  5-45.
Note      : converted here from mol/kg-H2O to mol/kg-SW
"""
@inline function Fᵈⁱᶜₖ₁ᵣ₉₃(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, p₁ = Pᵈⁱᶜₖ₁ᵣ₉₃, p₂ = Pᴴ²⁰ˢʷ) :: Real
    (; a₀, a₁, a₂, b₀, b₁, b₂, b₃, v₀, v₁, v₂, k₀, k₁, k₂) = p₁
    ln_kc1_p0 = (a₀ + 
                 a₁/Θᴷ + 
                 a₂*log(Θᴷ)
                + (
                   b₀ + 
                   b₁/Θᴷ
                  )*Sᴾ⁰⁵(Sᵖ)
                + b₂*Sᵖ
                + b₃*Sᵖ*Sᴾ⁰⁵(Sᵖ)
               )

    "Pressure correction for applied pressure /= 0"
    ln_kc1_pp = (-(v₀ + 
                   v₁*ΔSᵖ(Sᵖ) +
                   v₂*ΘᵒC(Θᴷ)) 
                 +(
                    k₀ + 
                    k₁*ΔSᵖ(Sᵖ) +
                    k₂*ΘᵒC(Θᴷ)
                )*(Δpᵦₐᵣ/2.)
                )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    return exp(ln_kc1_p0 + ln_kc1_pp) * H₂Oˢʷ(Sᵖ, p₂)
end

# Base.@kwdef struct Pᵈⁱᶜₖ₂ᵣ₉₃{FT}
#     a₀ :: FT = -   9.226508
#     a₁ :: FT = -3351.6106
#     a₂ :: FT = -   0.2005743
#     b₀ :: FT = -   0.106901773
#     b₁ :: FT = -  23.9722
#     b₂ :: FT =     0.1130822
#     b₃ :: FT = -   0.00846934
#     v₀ :: FT = -  15.82
#     v₁ :: FT =     0.321
#     v₂ :: FT = -   0.0219
#     k₀ :: FT =     1.13e-3
#     k₁ :: FT = -   0.314e-3
#     k₂ :: FT = -   0.1475e-3
# end
# adapt_structure( 
#         to, p::Pᵈⁱᶜₖ₂ᵣ₉₃,
#     ) = Pᵈⁱᶜₖ₂ᵣ₉₃(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
#     adapt(to, p.b₀),
#     adapt(to, p.b₁),
#     adapt(to, p.b₂),
#     adapt(to, p.b₃),
#     adapt(to, p.v₀),
#     adapt(to, p.v₁),
#     adapt(to, p.v₂),
#     adapt(to, p.k₀),
#     adapt(to, p.k₁),
#     adapt(to, p.k₂),
# )
"""
    Fᵈⁱᶜₖ₂ᵣ₉₃(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, Pᵈⁱᶜₖ₂ᵣ₉₃)

Return the second dissociation constant of carbonic acid in mol/kg-SW, given
temperature in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and
coefficients, `Pᵈⁱᶜₖ₂ᵣ₉₃`.

HCO₃⁻ <-> CO₃²⁻ + H⁺

References: Roy et al. (1993) -- also Handbook (1994)
            Millero (1979) pressure correction
pH scale  : Total
Valid range: T:  0-45  S:  5-45.
Note      : converted here from mol/kg-H2O to mol/kg-SW
"""
@inline function Fᵈⁱᶜₖ₂ᵣ₉₃(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, p₁ = Pᵈⁱᶜₖ₂ᵣ₉₃, p₂ = Pᴴ²⁰ˢʷ) :: Real
    (; a₀, a₁, a₂, b₀, b₁, b₂, b₃, v₀, v₁, v₂, k₀, k₁, k₂) = p₁
    ln_kc2_p0 = (a₀ + 
                 a₁/Θᴷ + 
                 a₂*log(Θᴷ)
                + (
                   b₀ + 
                   b₁/Θᴷ
                  )*Sᴾ⁰⁵(Sᵖ)
                + b₂*Sᵖ
                + b₃*Sᵖ*Sᴾ⁰⁵(Sᵖ)
               )    

    "Pressure correction for applied pressure /= 0"
    ln_kc2_pp = (-(v₀ + 
                   v₁*ΔSᵖ(Sᵖ) +
                   v₂*ΘᵒC(Θᴷ)) 
                 +(
                   k₀ + 
                   k₁*ΔSᵖ(Sᵖ) +
                   k₂*ΘᵒC(Θᴷ)
                )*(Δpᵦₐᵣ/2.)
                )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    return exp(ln_kc2_p0 + ln_kc2_pp) * H₂Oˢʷ(Sᵖ, p₂)
end

# Base.@kwdef struct Pᵈⁱᶜₖ₁ₘ₉₅{FT}
#     a₀ :: FT =     2.18867
#     a₁ :: FT = -2275.0360
#     a₂ :: FT = -   1.468591
#     b₀ :: FT = -   0.138681
#     b₁ :: FT = -   9.33291
#     b₂ :: FT =     0.0726483
#     b₃ :: FT = -   0.00574938
#     v₀ :: FT = -  25.5
#     v₁ :: FT = -   0.151
#     v₂ :: FT =     0.1271
#     k₀ :: FT = -   3.08e-3
#     k₁ :: FT = -   0.578e-3
#     k₂ :: FT =     0.0877e-3
# end
# adapt_structure( 
#         to, p::Pᵈⁱᶜₖ₁ₘ₉₅,
#     ) = Pᵈⁱᶜₖ₁ₘ₉₅(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
#     adapt(to, p.b₀),
#     adapt(to, p.b₁),
#     adapt(to, p.b₂),
#     adapt(to, p.b₃),
#     adapt(to, p.v₀),
#     adapt(to, p.v₁),
#     adapt(to, p.v₂),
#     adapt(to, p.k₀),
#     adapt(to, p.k₁),
#     adapt(to, p.k₂),
# )
"""
    Fᵈⁱᶜₖ₁ₘ₉₅(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, Pᵈⁱᶜₖ₁ₘ₉₅)

Return the first dissociation constant of carbonic acid in mol/kg-SW, given
temperature in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, 
`Δpᵦₐᵣ`, and coefficients, `Pᵈⁱᶜₖ₁ₘ₉₅`.
    
H₂CO₃ <-> HCO₃⁻ + H⁺
    
References: Millero (1995, eq 50 -- ln K1(COM))
             Millero (1982) pressure correction
pH scale:   SWS
"""
@inline function Fᵈⁱᶜₖ₁ₘ₉₅(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, p₁ = Pᵈⁱᶜₖ₁ₘ₉₅) :: Real
    (; a₀, a₁, a₂, b₀, b₁, b₂, b₃, v₀, v₁, v₂, k₀, k₁, k₂) = p₁
    ln_kc1_p0 = (a₀ + 
                 a₁/Θᴷ + 
                 a₂*log(Θᴷ)
                + (
                   b₀ + 
                   b₁/Θᴷ
                  )*Sᴾ⁰⁵(Sᵖ)
                + b₂*Sᵖ
                + b₃*Sᵖ*Sᴾ⁰⁵(Sᵖ)
               )


    "Pressure correction for applied pressure /= 0"
    ln_kc1_pp = (-(v₀ + 
                   v₁*ΔSᵖ(Sᵖ) +
                   v₂*ΘᵒC(Θᴷ)) 
                 +(
                   k₀ + 
                   k₁*ΔSᵖ(Sᵖ) +
                   k₂*ΘᵒC(Θᴷ)
                )*(Δpᵦₐᵣ/2.)
                )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    return exp(ln_kc1_p0 + ln_kc1_pp)
end

# Base.@kwdef struct Pᵈⁱᶜₖ₂ₘ₉₅{FT}
#     a₀ :: FT = -   0.84226
#     a₁ :: FT = -3741.1288
#     a₂ :: FT = -   1.437139
#     b₀ :: FT = -   0.128417
#     b₁ :: FT = -  24.41239
#     b₂ :: FT =     0.1195308
#     b₃ :: FT = -   0.00912840
#     v₀ :: FT = -  15.82
#     v₁ :: FT =     0.321
#     v₂ :: FT = -   0.0219
#     k₀ :: FT =     1.13e-3
#     k₁ :: FT = -   0.314e-3
#     k₂ :: FT = -   0.1475e-3
# end
# adapt_structure( 
#         to, p::Pᵈⁱᶜₖ₂ₘ₉₅,
#     ) = Pᵈⁱᶜₖ₂ₘ₉₅(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
#     adapt(to, p.b₀),
#     adapt(to, p.b₁),
#     adapt(to, p.b₂),
#     adapt(to, p.b₃),
#     adapt(to, p.v₀),
#     adapt(to, p.v₁),
#     adapt(to, p.v₂),
#     adapt(to, p.k₀),
#     adapt(to, p.k₁),
#     adapt(to, p.k₂),
# )
"""
    Fᵈⁱᶜₖ₂ₘ₉₅(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, Pᵈⁱᶜₖ₂ₘ₉₅)

Return the second dissociation constant of carbonic acid in mol/kg-SW, given
temperature in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and
coefficients, `Pᵈⁱᶜₖ₂ₘ₉₅`.
    
HCO₃⁻ <-> CO₃²⁻ + H⁺

References: Millero (1995, eq 51 -- ln K2(COM))
            Millero (1979) pressure correction
pH scale:   SWS
"""
@inline function Fᵈⁱᶜₖ₂ₘ₉₅(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, p₁ = Pᵈⁱᶜₖ₂ₘ₉₅) :: Real
    (; a₀, a₁, a₂, b₀, b₁, b₂, b₃, v₀, v₁, v₂, k₀, k₁, k₂) = p₁
    ln_kc2_p0 = (a₀ + 
                 a₁/Θᴷ + 
                 a₂*log(Θᴷ)
                + (
                   b₀ + 
                   b₁/Θᴷ
                )*Sᴾ⁰⁵(Sᵖ) +
                 b₂*Sᵖ +
                 b₃*Sᵖ*Sᴾ⁰⁵(Sᵖ)
                )

    "Pressure correction for applied pressure /= 0"
    ln_kc2_pp = (-(v₀ + 
                   v₁*ΔSᵖ(Sᵖ) +
                   v₂*ΘᵒC(Θᴷ)) 
                 +(
                   k₀ + 
                   k₁*ΔSᵖ(Sᵖ) +
                   k₂*ΘᵒC(Θᴷ)
                )*(Δpᵦₐᵣ/2.)
            )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    return exp(ln_kc2_p0 + ln_kc2_pp)
end

# Base.@kwdef struct Pᵈⁱᶜₖ₁ₗ₀₀{FT}
#     a₀ :: FT =    61.2172
#     a₁ :: FT = -3633.86
#     a₂ :: FT = -   9.67770
#     b₀ :: FT =     0.011555
#     b₁ :: FT = -   0.0001152
#     v₀ :: FT = -  25.5
#     v₁ :: FT = -   0.151
#     v₂ :: FT =     0.1271
#     k₀ :: FT = -   3.08e-3
#     k₁ :: FT = -   0.578e-3
#     k₂ :: FT =     0.0877e-3
# end
# adapt_structure( 
#         to, p::Pᵈⁱᶜₖ₁ₗ₀₀,
#     ) = Pᵈⁱᶜₖ₁ₗ₀₀(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
#     adapt(to, p.b₀),
#     adapt(to, p.b₁),
#     adapt(to, p.v₀),
#     adapt(to, p.v₁),
#     adapt(to, p.v₂),
#     adapt(to, p.k₀),
#     adapt(to, p.k₁),
#     adapt(to, p.k₂),
# )
"""
    Fᵈⁱᶜₖ₁ₗ₀₀(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, Pᵈⁱᶜₖ₁ₗ₀₀)

Return the first dissociation constant of carbonic acid in mol/kg-SW, given
temperature in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, 
`Δpᵦₐᵣ`, and coefficients, `Pᵈⁱᶜₖ₁ₗ₀₀`.
    
H₂CO₃ <-> HCO₃⁻ + H⁺
    
References: Luecker et al. (2000) -- also Handbook (2007)
            Millero (1979) pressure correction
pH scale:   Total
"""
@inline function Fᵈⁱᶜₖ₁ₗ₀₀(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, p₁ = Pᵈⁱᶜₖ₁ₗ₀₀) :: Real
    (; a₀, a₁, a₂, b₀, b₁, v₀, v₁, v₂, k₀, k₁, k₂) = p₁
    log10_kc1_p0 = (a₀ + 
                    a₁/Θᴷ + 
                    a₂*log(Θᴷ)
                   + (
                      b₀ + 
                      b₁*Sᵖ
                   )*Sᵖ
                  )

    "Pressure correction for applied pressure /= 0"
    ln_kc1_pp = (-(v₀ + 
                   v₁*ΔSᵖ(Sᵖ) +
                   v₂*ΘᵒC(Θᴷ)) 
                 +(
                   k₀ + 
                   k₁*ΔSᵖ(Sᵖ) +
                   k₂*ΘᵒC(Θᴷ)
                )*(Δpᵦₐᵣ/2.)
            )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    return 10^(log10_kc1_p0) * exp(ln_kc1_pp)
end

# Base.@kwdef struct Pᵈⁱᶜₖ₂ₗ₀₀{FT}
#     a₀ :: FT = -  25.9290
#     a₁ :: FT = - 471.78
#     a₂ :: FT =     3.16967
#     b₀ :: FT =     0.01781
#     b₁ :: FT = -   0.0001122
#     v₀ :: FT = -  15.82
#     v₁ :: FT =     0.321
#     v₂ :: FT = -   0.0219
#     k₀ :: FT =     1.13e-3
#     k₁ :: FT = -   0.314e-3
#     k₂ :: FT = -   0.1475e-3
# end
# adapt_structure( 
#         to, p::Pᵈⁱᶜₖ₂ₗ₀₀,
#     ) = Pᵈⁱᶜₖ₂ₗ₀₀(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
#     adapt(to, p.b₀),
#     adapt(to, p.b₁),
#     adapt(to, p.v₀),
#     adapt(to, p.v₁),
#     adapt(to, p.v₂),
#     adapt(to, p.k₀),
#     adapt(to, p.k₁),
#     adapt(to, p.k₂),
# )
"""
    Fᵈⁱᶜₖ₂ₗ₀₀(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, Pᵈⁱᶜₖ₂ₗ₀₀)

Return the second dissociation constant of carbonic acid in mol/kg-SW, given
temperature in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and
coefficients, `Pᵈⁱᶜₖ₂ₗ₀₀`.
    
HCO₃⁻ <-> CO₃²⁻ + H⁺

References: Luecker et al. (2000) -- also Handbook (2007)
            Millero (1979) pressure correction
pH scale:   Total
"""
@inline function Fᵈⁱᶜₖ₂ₗ₀₀(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, p₁ = Pᵈⁱᶜₖ₂ₗ₀₀) :: Real
    (; a₀, a₁, a₂, b₀, b₁, v₀, v₁, v₂, k₀, k₁, k₂) = p₁
    log10_kc2_p0 = (a₀ + 
                    a₁/Θᴷ + 
                    a₂*log(Θᴷ)
                   + (
                      b₀ + 
                      b₁*Sᵖ
                   )*Sᵖ
                  )

    "Pressure correction for applied pressure /= 0"
    ln_kc2_pp = (-(v₀ + 
                   v₁*ΔSᵖ(Sᵖ) +
                   v₂*ΘᵒC(Θᴷ)) 
                 +(
                   k₀ + 
                   k₁*ΔSᵖ(Sᵖ) +
                   k₂*ΘᵒC(Θᴷ)
                )*(Δpᵦₐᵣ/2.)
            )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    return 10^(log10_kc2_p0) * exp(ln_kc2_pp)
end

# Base.@kwdef struct Pᴮₖ₁{FT}
#     a₀ :: FT = - 8966.90
#     a₁ :: FT = - 2890.53
#     a₂ :: FT = -   77.942
#     a₃ :: FT =      1.728
#     a₄ :: FT = -    0.0996
#     b₀ :: FT =    148.0248
#     b₁ :: FT =    137.1942
#     b₂ :: FT =      1.62142
#     c₀ :: FT = -   24.4344
#     c₁ :: FT = -   25.085
#     c₂ :: FT = -    0.2474
#     d₀ :: FT =      0.053105
#     v₀ :: FT = -   29.48
#     v₁ :: FT =      0.295
#     v₂ :: FT =      0.1622
#     v₃ :: FT = -    0.002608
#     k₀ :: FT = -    2.84e-3
#     k₁ :: FT =      0.354e-3
# end
# adapt_structure( 
#         to, p::Pᴮₖ₁,
#     ) = Pᴮₖ₁(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
#     adapt(to, p.a₃),
#     adapt(to, p.a₄),
#     adapt(to, p.b₀),
#     adapt(to, p.b₁),
#     adapt(to, p.b₂),
#     adapt(to, p.c₀),
#     adapt(to, p.c₁),
#     adapt(to, p.c₂),
#     adapt(to, p.d₀),
#     adapt(to, p.v₀),
#     adapt(to, p.v₁),
#     adapt(to, p.v₂),
#     adapt(to, p.v₃),
#     adapt(to, p.k₀),
#     adapt(to, p.k₁),
# )
"""
    Fᵇₖ₁(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, Pᴮₖ₁)

Return boric acid dissociation constant in mol/kg-SW, given temperature in K,
`Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᴮₖ₁`.

References: Dickson (1990, eq. 23) -- also Handbook (2007, eq. 37)
            Millero (1979) pressure correction
pH scale  : total
"""
@inline function Fᵇₖ₁(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, p₁ = Pᴮₖ₁) :: Real
    (; a₀, a₁, a₂, a₃, a₄, b₀, b₁, b₂, c₀, c₁, c₂, d₀, v₀, v₁, v₂, v₃, k₀, k₁) = p₁
    ln_kb_p0  = ((a₀ +
                Sᴾ⁰⁵(Sᵖ)*(a₁ +
                Sᴾ⁰⁵(Sᵖ)*(a₂ +
                Sᴾ⁰⁵(Sᵖ)*(a₃ + 
                Sᴾ⁰⁵(Sᵖ)*  a₄
                   )))) / Θᴷ
                + b₀ + 
                Sᴾ⁰⁵(Sᵖ)*(
                    b₁ + 
                    b₂*Sᴾ⁰⁵(Sᵖ)
               ) +
                (c₀ + 
                Sᴾ⁰⁵(Sᵖ)*(
                    c₁ +
                    c₂*Sᴾ⁰⁵(Sᵖ))
                   ) * log(Θᴷ)
                + d₀*Sᴾ⁰⁵(Sᵖ)*Θᴷ
               )

    "Pressure correction for applied pressure /= 0"
    ln_kb_pp = (-(v₀ + 
                  v₁*ΔSᵖ(Sᵖ) +
                  v₂*ΘᵒC(Θᴷ) + 
                  v₃*ΘᵒC(Θᴷ)*ΘᵒC(Θᴷ)) 
                 +(
                   k₀ + 
                   k₁*ΔSᵖ(Sᵖ)
                )*(Δpᵦₐᵣ/2.)
                )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    return exp(ln_kb_p0 + ln_kb_pp)
end

# Base.@kwdef struct Pᴴ²ᴼₖ₁{FT}
#     a₀ :: FT =    148.9802
#     a₁ :: FT = -13847.26
#     a₂ :: FT = -   23.6521
#     b₀ :: FT = -    5.977
#     b₁ :: FT =    118.67
#     b₂ :: FT =      1.0495
#     c₀ :: FT = -    0.01615
#     v₀ :: FT =  -  20.02
#     v₁ :: FT =      0.1119
#     v₂ :: FT =  -   0.1409E-02
#     k₀ :: FT = -    5.13e-3
#     k₁ :: FT =      0.0794e-3
# end
# adapt_structure( 
#         to, p::Pᴴ²ᴼₖ₁,
#     ) = Pᴴ²ᴼₖ₁(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
#     adapt(to, p.b₀),
#     adapt(to, p.b₁),
#     adapt(to, p.b₂),
#     adapt(to, p.c₀),
#     adapt(to, p.v₀),
#     adapt(to, p.v₁),
#     adapt(to, p.v₂),
#     adapt(to, p.k₀),
#     adapt(to, p.k₁),
# )
"""
    Fᴴ²ᴼₖ₁(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, Pᴴ²ᴼₖ₁)

Return dissociation constant of water in (mol/kg-SW)^2, given temperature in K,
`Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᴴ²ᴼₖ₁`.

References: Millero (1995) for value at p_bar = 0
            Millero (pers. comm. 1996) for pressure correction
pH scale  : SWS
"""
@inline function Fᴴ²ᴼₖ₁(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, p₁ = Pᴴ²ᴼₖ₁) :: Real
    (; a₀, a₁, a₂, b₀, b₁, b₂, c₀, v₀, v₁, v₂, k₀, k₁) = p₁
    ln_kw_p0 = (a₀ +
                a₁/Θᴷ +
                a₂*log(Θᴷ) +
                (b₀ +
                 b₁/Θᴷ +
                 b₂*log(Θᴷ)
               )*Sᴾ⁰⁵(Sᵖ) +
                 c₀*Sᵖ
               )

    "Pressure correction for applied pressure /= 0"
    ln_kw_pp = (-(v₀ + 
                  v₁*ΘᵒC(Θᴷ) + 
                  v₂*ΘᵒC(Θᴷ)*ΘᵒC(Θᴷ)) 
                +(
                  k₀ + 
                  k₁*ΘᵒC(Θᴷ)
                )*(Δpᵦₐᵣ/2.)
                )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    return exp(ln_kw_p0 + ln_kw_pp)
end

# Base.@kwdef struct Pᴾᴼ⁴ₖ₁{FT}
#     a₀ :: FT =   115.54
#     a₁ :: FT = -4576.752
#     a₂ :: FT = -  18.453
#     b₀ :: FT =     0.69171
#     b₁ :: FT = - 106.736
#     b₂ :: FT = -   0.01844
#     b₃ :: FT = -   0.65643
#     v₀ :: FT = -  14.51
#     v₁ :: FT =     0.1211
#     v₂ :: FT = -   0.321E-03
#     k₀ :: FT = -   2.67e-3
#     k₁ :: FT =     0.0427e-3
# end
# adapt_structure( 
#         to, p::Pᴾᴼ⁴ₖ₁,
#     ) = Pᴾᴼ⁴ₖ₁(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
#     adapt(to, p.b₀),
#     adapt(to, p.b₁),
#     adapt(to, p.b₂),
#     adapt(to, p.b₃),
#     adapt(to, p.v₀),
#     adapt(to, p.v₁),
#     adapt(to, p.v₂),
#     adapt(to, p.k₀),
#     adapt(to, p.k₁),
# )
"""
    Fᴾᴼ⁴ₖ₁(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, Pᴾᴼ⁴ₖ₁)

Return the first dissociation constant of phosphoric acid (H3PO4) in seawater, given 
temperature in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᴾᴼ⁴ₖ₁`.

References: Yao and Millero (1995)
            Millero (1995) for pressure correction
pH scale  : SWS
"""
@inline function Fᴾᴼ⁴ₖ₁(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, p₁ = Pᴾᴼ⁴ₖ₁) :: Real
    (; a₀, a₁, a₂, b₀, b₁, b₂, b₃, v₀, v₁, v₂, k₀, k₁) = p₁
    ln_kp1_p0 = (a₀ +
                 a₁/Θᴷ +
                 a₂*log(Θᴷ) +
                (
                 b₀ +
                 b₁/Θᴷ
                )*Sᴾ⁰⁵(Sᵖ) +
                 (
                  b₂ +
                  b₃/Θᴷ
                )*Sᵖ)

    "Pressure correction for applied pressure /= 0"
    ln_kp1_pp = (-(v₀ + 
                   v₁*ΘᵒC(Θᴷ) +
                   v₂*ΘᵒC(Θᴷ)*ΘᵒC(Θᴷ)) 
                 +(
                   k₀ + 
                   k₁*ΘᵒC(Θᴷ)
                )*(Δpᵦₐᵣ/2.)
                )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    return exp(ln_kp1_p0 + ln_kp1_pp)
end

# Base.@kwdef struct Pᴾᴼ⁴ₖ₂{FT}
#     a₀ :: FT =    172.1033
#     a₁ :: FT = - 8814.715
#     a₂ :: FT = -   27.927
#     b₀ :: FT =      1.3566
#     b₁ :: FT = -  160.340
#     b₂ :: FT = -    0.05778
#     b₃ :: FT =      0.37335
#     v₀ :: FT = -   23.12
#     v₁ :: FT =      0.1758
#     v₂ :: FT = -    0.002647
#     k₀ :: FT = -    5.15e-3
#     k₁ :: FT =      0.09e-3
# end
# adapt_structure( 
#         to, p::Pᴾᴼ⁴ₖ₂,
#     ) = Pᴾᴼ⁴ₖ₂(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
#     adapt(to, p.b₀),
#     adapt(to, p.b₁),
#     adapt(to, p.b₂),
#     adapt(to, p.b₃),
#     adapt(to, p.v₀),
#     adapt(to, p.v₁),
#     adapt(to, p.v₂),
#     adapt(to, p.k₀),
#     adapt(to, p.k₁),
# )
"""
    Fᴾᴼ⁴ₖ₂(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, Pᴾᴼ⁴ₖ₂)

Return the second dissociation constant of phosphoric acid (H3PO4) in seawater, given
temperature in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᴾᴼ⁴ₖ₂`.

References: Yao and Millero (1995)
            Millero (1995) for pressure correction
pH scale  : SWS
"""
@inline function Fᴾᴼ⁴ₖ₂(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, p₁ = Pᴾᴼ⁴ₖ₂) :: Real
    (; a₀, a₁, a₂, b₀, b₁, b₂, b₃, v₀, v₁, v₂, k₀, k₁) = p₁
    ln_kp2_p0 = (a₀ +
                 a₁/Θᴷ +
                 a₂*log(Θᴷ) +
                (
                 b₀ +
                 b₁/Θᴷ
                )*Sᴾ⁰⁵(Sᵖ) +
                 (
                  b₂ +
                  b₃/Θᴷ
                )*Sᵖ)

    "Pressure correction for applied pressure /= 0"
    ln_kp2_pp = (-(v₀ + 
                   v₁*ΘᵒC(Θᴷ) +
                   v₂*ΘᵒC(Θᴷ)*ΘᵒC(Θᴷ)) 
                 +(
                   k₀ + 
                   k₁*ΘᵒC(Θᴷ)
                )*(Δpᵦₐᵣ/2.)
                )*(Δpᵦₐᵣ/Rₜ(Θᴷ))
    
    return exp(ln_kp2_p0 + ln_kp2_pp)
end

# Base.@kwdef struct Pᴾᴼ⁴ₖ₃{FT}
#     a₀ :: FT = -   18.126
#     a₁ :: FT = - 3070.75
#     a₂ :: FT =      2.81197
#     a₃ :: FT =     17.27039
#     a₄ :: FT = -    0.09984
#     a₅ :: FT = -   44.99486
#     v₀ :: FT = -   26.57
#     v₁ :: FT =      0.2020
#     v₂ :: FT = -    3.042e-3
#     k₀ :: FT = -    4.08e-3
#     k₁ :: FT =      0.0714e-3
# end
# adapt_structure( 
#         to, p::Pᴾᴼ⁴ₖ₃,
#     ) = Pᴾᴼ⁴ₖ₃(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
#     adapt(to, p.a₃),
#     adapt(to, p.a₄),
#     adapt(to, p.a₅),
#     adapt(to, p.v₀),
#     adapt(to, p.v₁),
#     adapt(to, p.v₂),
#     adapt(to, p.k₀),
#     adapt(to, p.k₁),
# )
"""
    Fᴾᴼ⁴ₖ₃(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, Pᴾᴼ⁴ₖ₃)

Return the third dissociation constant of phosphoric acid (H3PO4) in seawater, given 
temperature in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᴾᴼ⁴ₖ₃`.

References: Yao and Millero (1995)
            Millero (1995) for pressure correction
pH scale  : SWS
"""
@inline function Fᴾᴼ⁴ₖ₃(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, p₁ = Pᴾᴼ⁴ₖ₃) :: Real
    (; a₀, a₁, a₂, a₃, a₄, a₅, v₀, v₁, v₂, k₀, k₁) = p₁
    ln_kp3_p0 = (a₀ +
                 a₁/Θᴷ +
                (
                 a₂ +
                 a₃/Θᴷ
                )*Sᴾ⁰⁵(Sᵖ) +
                (
                 a₄ +
                 a₅/Θᴷ
                )*Sᵖ)

    "Pressure correction for applied pressure /= 0"
    ln_kp3_pp = (-(v₀ + 
                   v₁*ΘᵒC(Θᴷ) +
                   v₂*ΘᵒC(Θᴷ)*ΘᵒC(Θᴷ)) 
                 +(
                   k₀ + 
                   k₁*ΘᵒC(Θᴷ)
                )*(Δpᵦₐᵣ/2.)
                )*(Δpᵦₐᵣ/Rₜ(Θᴷ))
    
    return exp(ln_kp3_p0 + ln_kp3_pp)
end

# Base.@kwdef struct Pˢⁱᵗₖ₁{FT}
#     a₀ :: FT =    117.40
#     a₁ :: FT = - 8904.2
#     a₂ :: FT = -   19.334
#     b₀ :: FT =      3.5913
#     b₁ :: FT = -  458.79
#     b₂ :: FT = -    1.5998
#     b₃ :: FT =    188.74
#     c₀ :: FT =      0.07871
#     c₁ :: FT = -   12.1652
#     v₀ :: FT = -   29.48
#     v₁ :: FT =      0.0
#     v₂ :: FT =      0.1622
#     v₃ :: FT = -    0.002608
#     k₀ :: FT = -    2.84e-3
#     k₁ :: FT =      0.354e-3
# end
# adapt_structure( 
#         to, p::Pˢⁱᵗₖ₁,
#     ) = Pˢⁱᵗₖ₁(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
#     adapt(to, p.b₀),
#     adapt(to, p.b₁),
#     adapt(to, p.b₂),
#     adapt(to, p.b₃),
#     adapt(to, p.c₀),
#     adapt(to, p.c₁),
#     adapt(to, p.v₀),
#     adapt(to, p.v₁),
#     adapt(to, p.v₂),
#     adapt(to, p.v₃),
#     adapt(to, p.k₀),
#     adapt(to, p.k₁),
# )
"""
    Fˢⁱᵗₖ₁(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, Pˢⁱᵗₖ₁)

Return the first dissociation constant of silicic acid (H4SiO4) in seawater, given 
temperature in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pˢⁱᵗₖ₁`.

References: Yao and Millero (1995) cited by Millero (1995)
pH scale  : SWS (according to Dickson et al, 2007)
Note      : No pressure correction available
Note      : converted here from mol/kg-H2O to mol/kg-sw
"""
@inline function Fˢⁱᵗₖ₁(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, p₁ = Pˢⁱᵗₖ₁, p₂ = Pᵘˢ, p₃ = Pᴴ²⁰ˢʷ) :: Real
    (; a₀, a₁, a₂, b₀, b₁, b₂, b₃, c₀, c₁, v₀, v₁, v₂, v₃, k₀, k₁ ) = p₁

    ln_ksi1_p0 = (a₀ +
                  a₁/Θᴷ +
                  a₂*log(Θᴷ) +
                 (b₀ +
                  b₁/Θᴷ
                 )*sqrt(μₛ(Sᵖ, p₂, p₃)) +
                 (b₂ +
                  b₃/Θᴷ
                 )*μₛ(Sᵖ, p₂, p₃) +
                 (c₀ +
                  c₁/Θᴷ
                 )*μₛ(Sᵖ, p₂, p₃)*μₛ(Sᵖ, p₂, p₃))

    "Pressure correction : currently none"
    "ln_ksi1_pp = 0."
    "Pressure correction for applied pressure /= 0 
    estimated from borate as suggested by Millero (1995)"
    ln_kb_pp = (-(v₀ + 
                  v₁*ΔSᵖ(Sᵖ) +
                  v₂*ΘᵒC(Θᴷ) + 
                  v₃*ΘᵒC(Θᴷ)*ΘᵒC(Θᴷ)) 
                +(
                  k₀ + 
                  k₁*ΔSᵖ(Sᵖ)
                )*(Δpᵦₐᵣ/2.)
                )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    return exp(ln_ksi1_p0 + ln_kb_pp) * H₂Oˢʷ(Sᵖ, p₃)
end

# Base.@kwdef struct Pᴴ²ˢₖ₁{FT}
#     a₀ :: FT =     225.838
#     a₁ :: FT = - 13275.3
#     a₂ :: FT = -    34.6435
#     a₃ :: FT =       0.3449
#     a₄ :: FT = -     0.0274
#     v₀ :: FT = -    14.80
#     v₁ :: FT =       0.0020
#     v₂ :: FT =  -    0.400E-03
#     k₀ :: FT =       2.89e-3
#     k₁ :: FT =       0.054e-3
# end
# adapt_structure( 
#         to, p::Pᴴ²ˢₖ₁,
#     ) = Pᴴ²ˢₖ₁(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
#     adapt(to, p.a₃),
#     adapt(to, p.a₄),
#     adapt(to, p.v₀),
#     adapt(to, p.v₁),
#     adapt(to, p.v₂),
#     adapt(to, p.k₀),
#     adapt(to, p.k₁),
# )
"""
    Fᴴ²ˢₖ₁(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, Pᴴ²ˢₖ₁)

Return the dissociation constant of hydrogen sulfide in sea-water, given temperature in K,
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
@inline function Fᴴ²ˢₖ₁(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, p₁ = Pᴴ²ˢₖ₁) :: Real
    (; a₀, a₁, a₂, a₃, a₄, v₀, v₁, v₂, k₀, k₁) = p₁
    ln_kh2s_p0 = (a₀ +
                  a₁/Θᴷ +
                  a₂*log(Θᴷ) +
                  a₃*Sᴾ⁰⁵(Sᵖ) +
                  a₄*Sᵖ)

    "Pressure correction for applied pressure /= 0"
    ln_kh2s_pp = (-(v₀ +
                    v₁*ΘᵒC(Θᴷ) +
                    v₂*ΘᵒC(Θᴷ)*ΘᵒC(Θᴷ)
                  )+(
                     k₀ +
                     k₁*ΘᵒC(Θᴷ)
                    )*(Δpᵦₐᵣ/2.)
                    )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    return exp(ln_kh2s_p0 + ln_kh2s_pp)
end

# Base.@kwdef struct Pᴺᴴ⁴ₖ₁{FT}
#     a₀ :: FT = -    0.25444
#     a₁ :: FT = - 6285.33
#     a₂ :: FT =      0.0001635
#     b₀ :: FT =      0.46532
#     b₁ :: FT = -  123.7184
#     b₂ :: FT = -    0.01992
#     b₃ :: FT =      3.17556
#     v₀ :: FT = -   26.43
#     v₁ :: FT =      0.0889
#     v₂ :: FT = -    0.905E-03
#     k₀ :: FT = -    5.03E-03
#     k₁ :: FT =      0.0814E-03
# end
# adapt_structure( 
#         to, p::Pᴺᴴ⁴ₖ₁,
#     ) = Pᴺᴴ⁴ₖ₁(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
#     adapt(to, p.b₀),
#     adapt(to, p.b₁),
#     adapt(to, p.b₂),
#     adapt(to, p.b₃),
#     adapt(to, p.v₀),
#     adapt(to, p.v₁),
#     adapt(to, p.v₂),
#     adapt(to, p.k₀),
#     adapt(to, p.k₁),
# )
"""
    Fᴺᴴ⁴ₖ₁(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, Pᴺᴴ⁴ₖ₁)

Return the dissociation constant of ammonium in sea-water [mol/kg-SW], given temperature
in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᴺᴴ⁴ₖ₁`.

References: Yao and Millero (1995)
            Millero (1995) for pressure correction
pH scale  : SWS
"""
@inline function Fᴺᴴ⁴ₖ₁(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, p₁ = Pᴺᴴ⁴ₖ₁) :: Real
    (; a₀, a₁, a₂, b₀, b₁, b₂, b₃, v₀, v₁, v₂, k₀, k₁) = p₁
    ln_knh4_p0 = (a₀ +
                  a₁/Θᴷ +
                  a₂*Θᴷ +
                  (b₀ +
                   b₁/Θᴷ
                   )*Sᴾ⁰⁵(Sᵖ) +
                    (b₂ +
                     b₃/Θᴷ
                    )*Sᵖ)

    
    "Pressure correction for applied pressure /= 0"
    ln_knh4_pp = (-(v₀ +
                    v₁*ΘᵒC(Θᴷ) +
                    v₂*ΘᵒC(Θᴷ)*ΘᵒC(Θᴷ)
                  )+(
                    k₀ +
                    k₁*ΘᵒC(Θᴷ)
                    )*(Δpᵦₐᵣ/2.)
                    )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    return exp(ln_knh4_p0 + ln_knh4_pp)
end

# Base.@kwdef struct Pᴴᶠᵦ₁{FT}
#     a₀ :: FT =     12.641
#     a₁ :: FT = - 1590.2
#     a₂ :: FT = -    1.525
#     v₀ :: FT = -    9.78
#     v₁ :: FT = -    0.0090
#     v₂ :: FT = -    0.942E-03
#     k₀ :: FT = -    3.91e-3
#     k₁ :: FT =      0.054e-3
# end
# adapt_structure( 
#         to, p::Pᴴᶠᵦ₁,
#     ) = Pᴴᶠᵦ₁(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
#     adapt(to, p.v₀),
#     adapt(to, p.v₁),
#     adapt(to, p.v₂),
#     adapt(to, p.k₀),
#     adapt(to, p.k₁),
# )
"""
    Fᴴᶠᵦ₁(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, Pᴴᶠᵦ₁)

Return the association constant of HF in (mol/kg-SW)^-1, given temperature
in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᴴᶠᵦ₁`.

HF <-> H⁺ + F⁻ 

References: Dickson and Riley (1979)
            Millero (1995) for pressure correction
pH scale  : free
Note      : converted here from mol/kg-H2O to mol/kg-SW
"""
@inline function Fᴴᶠᵦ₁(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, p₁ = Pᴴᶠᵦ₁, p₂ = Pᵘˢ, p₃ = Pᴴ²⁰ˢʷ) :: Real
    (; a₀, a₁, a₂, v₀, v₁, v₂, k₀, k₁) = p₁
    ln_bhf_p0 = (a₀ +
                 a₁/Θᴷ +
                 a₂*sqrt(μₛ(Sᵖ, p₂, p₃)))

    "Pressure correction for applied pressure /= 0"
    ln_khf_pp = (-(v₀ +
                   v₁*ΘᵒC(Θᴷ) +
                   v₂*ΘᵒC(Θᴷ)*ΘᵒC(Θᴷ)
                  )+(
                   k₀ +
                   k₁*ΘᵒC(Θᴷ)
                   )*(Δpᵦₐᵣ/2.)
                   )*(Δpᵦₐᵣ/Rₜ(Θᴷ))
    """                
    Final \beta_HF value
    --------------------
     notice that  ln(k_HF(P)) = ln(k_HF(0)) + zln_khf_pp
            <=>  -ln(\beta_HF(P)) = -ln(\beta_HF(0)) + zln_khf_pp
            <=>   ln(\beta_HF(P)) =  ln(\beta_HF(0)) - zln_khf_pp
    """
    return exp(ln_bhf_p0 - ln_khf_pp) / H₂Oˢʷ(Sᵖ, p₃)
end

# Base.@kwdef struct Pᴴᶠₖ₁{FT}
#     a₀ :: FT = -  9.68
#     a₁ :: FT =  874.0
#     a₂ :: FT =    0.111
#     v₀ :: FT = -  9.78
#     v₁ :: FT = -  0.0090
#     v₂ :: FT = -  0.942E-3
#     k₀ :: FT = -  3.91e-3
#     k₁ :: FT =    0.054e-3
# end
# adapt_structure( 
#         to, p::Pᴴᶠₖ₁,
#     ) = Pᴴᶠₖ₁(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
#     adapt(to, p.v₀),
#     adapt(to, p.v₁),
#     adapt(to, p.v₂),
#     adapt(to, p.k₀),
#     adapt(to, p.k₁),
# )
"""
    Fᴴᶠₖ₁(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, Pᴴᶠₖ₁)

Return the dissociation constant for hydrogen fluoride in mol/kg-SW, given temperature
in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᴴᶠₖ₁`.

HF <-> H⁺ + F⁻ 

References: Perez and Fraga (1987)
            Millero (1995) for pressure correction
pH scale  : Total (according to Handbook, 2007
"""
@inline function Fᴴᶠₖ₁(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, p₁ = Pᴴᶠₖ₁) :: Real
    (; a₀, a₁, a₂, v₀, v₁, v₂, k₀, k₁) = p₁
    ln_khf_p0 = (a₀ +
                 a₁/Θᴷ +
                 a₂*Sᴾ⁰⁵(Sᵖ))

    "Pressure correction for applied pressure /= 0"
    ln_khf_pp = (-(v₀ +
                   v₁*ΘᵒC(Θᴷ) +
                   v₂*ΘᵒC(Θᴷ)*ΘᵒC(Θᴷ)
                  )+(
                    k₀ +
                    k₁*ΘᵒC(Θᴷ)
                    )*(Δpᵦₐᵣ/2.)
                    )*(Δpᵦₐᵣ/Rₜ(Θᴷ))
        
    return exp(ln_khf_p0 + ln_khf_pp)
end

# Base.@kwdef struct Pᴴˢᴼ⁴ₖ₁{FT}
#     a₀ :: FT =     141.328
#     a₁ :: FT = -  4276.1
#     a₂ :: FT = -    23.093
#     b₀ :: FT =     324.57
#     b₁ :: FT = - 13856.
#     b₂ :: FT = -    47.986
#     c₀ :: FT = -   771.54
#     c₁ :: FT =   35474.
#     c₂ :: FT =     114.723
#     d₀ :: FT = -  2698.
#     d₁ :: FT =    1776.
#     v₀ :: FT = -    18.03
#     v₁ :: FT =       0.0466
#     v₂ :: FT =       0.316E-03
#     k₀ :: FT = -     4.53e-3
#     k₁ :: FT =       0.0900e-3
# end
# adapt_structure( 
#         to, p::Pᴴˢᴼ⁴ₖ₁,
#     ) = Pᴴˢᴼ⁴ₖ₁(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
#     adapt(to, p.b₀),
#     adapt(to, p.b₁),
#     adapt(to, p.b₂),
#     adapt(to, p.c₀),
#     adapt(to, p.c₁),
#     adapt(to, p.c₂),
#     adapt(to, p.d₀),
#     adapt(to, p.d₁),
#     adapt(to, p.v₀),
#     adapt(to, p.v₁),
#     adapt(to, p.v₂),
#     adapt(to, p.k₀),
#     adapt(to, p.k₁),
# )
"""
    Fᴴˢᴼ⁴ₖ₁(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, Pᴴˢᴼ⁴ₖ₁)

Return the dissociation constant of hydrogen sulfate (bisulfate) , given temperature
in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᴴˢᴼ⁴ₖ₁`.

References: Dickson (1990) -- also Handbook (2007)
            Millero (1995) for pressure correction
pH scale  : free
Note      : converted here from mol/kg-H2O to mol/kg-SW
"""
@inline function Fᴴˢᴼ⁴ₖ₁(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, p₁ = Pᴴˢᴼ⁴ₖ₁, p₂ = Pᵘˢ, p₃ = Pᴴ²⁰ˢʷ) :: Real
    (; a₀, a₁, a₂, b₀, b₁, b₂, c₀, c₁, c₂, d₀, d₁, v₀, v₁, v₂, k₀, k₁) = p₁
    ln_khso4_p0 = (a₀ +
                   a₁/Θᴷ +
                   a₂*log(Θᴷ) +
                   (b₀ +
                    b₁/Θᴷ +
                    b₂*log(Θᴷ)
                  )*sqrt(μₛ(Sᵖ, p₂, p₃)) +
                   (c₀ +
                    c₁/Θᴷ +
                    c₂*log(Θᴷ)
                  )*μₛ(Sᵖ, p₂, p₃) +
                   (d₀/Θᴷ)*sqrt(μₛ(Sᵖ, p₂, p₃))*μₛ(Sᵖ, p₂, p₃)+
                   (d₁/Θᴷ)*μₛ(Sᵖ, p₂, p₃)*μₛ(Sᵖ, p₂, p₃))

    "Pressure correction for applied pressure /= 0"
    ln_khso4_pp = (-(v₀ +
                     v₁*ΘᵒC(Θᴷ) +
                     v₂*ΘᵒC(Θᴷ)*ΘᵒC(Θᴷ)
                   )+(
                      k₀ +
                      k₁*ΘᵒC(Θᴷ)
                     )*(Δpᵦₐᵣ/2.)
                     )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    return exp(ln_khso4_p0 + ln_khso4_pp) * H₂Oˢʷ(Sᵖ, p₃)
end

# Base.@kwdef struct Pᶜᵃˡᶜⁱᵗᵉₛₚ{FT}
#     a₀ :: FT = - 171.9065
#     a₁ :: FT = -   0.077993
#     a₂ :: FT =  2839.319
#     a₃ :: FT =    71.595
#     b₀ :: FT = -   0.77712
#     b₁ :: FT =     0.0028426
#     b₂ :: FT =   178.34
#     c₀ :: FT = -   0.07711
#     d₀ :: FT =     0.0041249
#     v₀ :: FT = -  48.76
#     v₁ :: FT =     0.5304
#     k₀ :: FT = -  11.76e-3
#     k₁ :: FT =     0.3692e-3
# end
# adapt_structure( 
#         to, p::Pᶜᵃˡᶜⁱᵗᵉₛₚ,
#     ) = Pᶜᵃˡᶜⁱᵗᵉₛₚ(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
#     adapt(to, p.a₃),
#     adapt(to, p.b₀),
#     adapt(to, p.b₁),
#     adapt(to, p.b₂),
#     adapt(to, p.c₀),
#     adapt(to, p.d₀),
#     adapt(to, p.v₀),
#     adapt(to, p.v₁),
#     adapt(to, p.k₀),
#     adapt(to, p.k₁),
# )
"""
    Fᶜᵃˡᶜⁱᵗᵉₛₚ(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, Pᶜᵃˡᶜⁱᵗᵉₛₚ)

Return the stoichiometric solubility product of calcite, `Ω`, in seawater, given temperature
in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᶜᵃˡᶜⁱᵗᵉₛₚ`

References: Mucci (1983)
            Millero (1995) for pressure correction
pH scale  : N/A
Units     : (mol/kg-SW)^2
"""
@inline function Fᶜᵃˡᶜⁱᵗᵉₛₚ(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, p₁ = Pᶜᵃˡᶜⁱᵗᵉₛₚ) :: Real
    (; a₀, a₁, a₂, a₃, b₀, b₁, b₂, c₀, d₀, v₀, v₁, k₀, k₁) = p₁
    log10_kcalcite_p0 = (a₀ +
                         a₁*Θᴷ +
                         a₂/Θᴷ +
                         a₃*log10(Θᴷ) +
                        (b₀ +
                         b₁*Θᴷ +
                         b₂/Θᴷ
                        )*Sᴾ⁰⁵(Sᵖ) +
                         (c₀)*Sᵖ +
                         (d₀)*Sᵖ*Sᴾ⁰⁵(Sᵖ))

    "Pressure correction for applied pressure /= 0"
    ln_kcalcite_pp = (-(v₀ +
                        v₁*ΘᵒC(Θᴷ)
                        )+(
                          k₀ +
                          k₁*ΘᵒC(Θᴷ)
                         )*(Δpᵦₐᵣ/2.)
                        )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    return 10^(log10_kcalcite_p0) * exp(ln_kcalcite_pp)
end

# Base.@kwdef struct Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ{FT}
#     a₀ :: FT = - 171.945
#     a₁ :: FT = -   0.077993
#     a₂ :: FT =  2903.293
#     a₃ :: FT =    71.595
#     b₀ :: FT = -   0.068393
#     b₁ :: FT =     0.0017276
#     b₂ :: FT =    88.135
#     c₀ :: FT = -   0.10018
#     d₀ :: FT =     0.0059415
#     v₀ :: FT = -  48.76
#     v₁ :: FT =     0.5304
#     v₂ :: FT =     2.8
#     k₀ :: FT = -  11.76e-3
#     k₁ :: FT =     0.3692e-3
# end
# adapt_structure( 
#         to, p::Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ,
#     ) = Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ(
#     adapt(to, p.a₀),
#     adapt(to, p.a₁),
#     adapt(to, p.a₂),
#     adapt(to, p.a₃),
#     adapt(to, p.b₀),
#     adapt(to, p.b₁),
#     adapt(to, p.b₂),
#     adapt(to, p.c₀),
#     adapt(to, p.d₀),
#     adapt(to, p.v₀),
#     adapt(to, p.v₁),
#     adapt(to, p.v₂),
#     adapt(to, p.k₀),
#     adapt(to, p.k₁),
# )
"""
    Fᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ)

Return stoichiometric solubility product, `Ω`, of aragonite in seawater, given temperature
in K, `Θᴷ`, practical salinity, `Sᵖ`, applied pressure, `Δpᵦₐᵣ`, and coefficients, `Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ`.

References: Mucci (1983)
            Millero (1979) for pressure correction
pH scale  : N/A
Units     : (mol/kg-SW)^2
"""
@inline function Fᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, p₁ = Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ) :: Real
    (; a₀, a₁, a₂, a₃, b₀, b₁, b₂, c₀, d₀, v₀, v₁, v₂, k₀, k₁) = p₁
    log10_karagonite_p0 = (a₀ +
                           a₁*Θᴷ +
                           a₂/Θᴷ +
                           a₃*log10(Θᴷ) +
                           (b₀ +
                            b₁*Θᴷ +
                            b₂/Θᴷ
                          )*Sᴾ⁰⁵(Sᵖ) +
                           (c₀)*Sᵖ +
                           (d₀)*Sᵖ*Sᴾ⁰⁵(Sᵖ))

    "Pressure correction for applied pressure /= 0"
    ln_karagonite_pp = (-(v₀ +
                          v₁*ΘᵒC(Θᴷ) +
                          v₂
                         )+(
                           k₀ +
                           k₁*ΘᵒC(Θᴷ)
                           )*(Δpᵦₐᵣ/2.)
                           )*(Δpᵦₐᵣ/Rₜ(Θᴷ))

    return 10^(log10_karagonite_p0) * exp(ln_karagonite_pp)
end

"""
    H⁺ₛoverH⁺ₜ(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real)

Return the ratio H_SWS/H_Tot as a function of salinity, `Sᵖ`.
Reference:  Munhoven
pH scale:   all
"""
@inline function H⁺ₛoverH⁺ₜ(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, p₁ = Pᶠᵀᴼᵀ, p₂ = Pᴴᶠᵦ₁, p₃ = Pᵘˢ, p₄ = Pᴴ²⁰ˢʷ, p₅ = Pˢᴼ⁴ᵀᴼᵀ, p₆ = Pᴴˢᴼ⁴ₖ₁) :: Real
    return (1. +  
            (Fᵀᴼᵀ(Sᵖ, p₁)*Fᴴᶠᵦ₁(Θᴷ, Sᵖ, Δpᵦₐᵣ, p₂, p₃, p₄))
           /(1. + SO₄ᵀᴼᵀ(Sᵖ, p₅)/Fᴴˢᴼ⁴ₖ₁(Θᴷ, Sᵖ, Δpᵦₐᵣ, p₆, p₃, p₄))
          )
end

"""
    H⁺ₜoverH⁺₃(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real)

Return the ratio H_Tot/H_free as a function of salinity, `Sᵖ`.
Reference:  Munhoven
pH scale:   N/A
"""
@inline function H⁺ₜoverH⁺₃(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, p₁ = Pˢᴼ⁴ᵀᴼᵀ, p₂ = Pᴴˢᴼ⁴ₖ₁, p₃ = Pᵘˢ, p₄ = Pᴴ²⁰ˢʷ) :: Real
    return 1. + SO₄ᵀᴼᵀ(Sᵖ, p₁)/Fᴴˢᴼ⁴ₖ₁(Θᴷ, Sᵖ, Δpᵦₐᵣ, p₂, p₃, p₄)
end

"""
    H⁺ₛoverH⁺₃(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real)

Return the ratio H_SWS/H_free as a function
of salinity, `Sᵖ`.
"""
@inline function H⁺ₛoverH⁺₃(Θᴷ::Real, Sᵖ::Real, Δpᵦₐᵣ::Real, p₁ = Pᶠᵀᴼᵀ, p₂ = Pᴴᶠᵦ₁, p₃ = Pᵘˢ, p₄ = Pᴴ²⁰ˢʷ, p₅ = Pˢᴼ⁴ᵀᴼᵀ, p₆ = Pᴴˢᴼ⁴ₖ₁) :: Real
    return (1. + 
            (Fᵀᴼᵀ(Sᵖ, p₁)*Fᴴᶠᵦ₁(Θᴷ, Sᵖ, Δpᵦₐᵣ, p₂, p₃, p₄)) +
            (SO₄ᵀᴼᵀ(Sᵖ, p₅)/Fᴴˢᴼ⁴ₖ₁(Θᴷ, Sᵖ, Δpᵦₐᵣ, p₆, p₃, p₄))
        )
end