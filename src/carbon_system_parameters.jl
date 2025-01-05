# Struct to hold parameters for the UniversalRobustCarbonSolver
Base.@kwdef struct CarbonSolverParameters #{FT<:Real, IT<:Int} #
    Δₕ₊      #:: FT # Tolerance of the H⁺ solution 
    H⁺ᵗʰʳᵉˢʰ #:: FT # H⁺ threshold for secant iteration
    Iᴴ⁺ₘₐₓ   #:: IT # Maximum number of iterations
end

"""
    CarbonSolverParameters(Δₕ₊::Real=1e-8, H⁺ᵗʰʳᵉˢʰ::Real=1, Iᴴ⁺ₘₐₓ::Real=100)

Create a `CarbonSolverParameters` object with the specified parameters.

# Arguments
- `grid`: The grid object which determines the element type for the parameters.
- `Δₕ₊`: A real number representing the increment for the hydrogen ion concentration. Default is `1e-8`.
- `H⁺ᵗʰʳᵉˢʰ`: A real number representing the threshold for the hydrogen ion concentration. Default is `1`.
- `Iᴴ⁺ₘₐₓ`: A real number representing the maximum number of iterations for the solver. Default is `100`.

# Returns
- A `CarbonSolverParameters` object with the specified parameters converted to the element type of the grid.
"""
function CarbonSolverParameters(;
    Δₕ₊      :: Real = 1e-8, 
    H⁺ᵗʰʳᵉˢʰ :: Real = 1.0, 
    Iᴴ⁺ₘₐₓ   :: Real = 100.0,
)
    return CarbonSolverParameters(
        Δₕ₊,
        H⁺ᵗʰʳᵉˢʰ,
        Iᴴ⁺ₘₐₓ,
    )
end

# Struct to hold carbon dissociation function coefficients for chemical species
Base.@kwdef struct CarbonCoefficientParameters #{FT<:Real}
    a₀ #:: FT
    a₁ #:: FT
    a₂ #:: FT
    a₃ #:: FT
    a₄ #:: FT
    a₅ #:: FT
    b₀ #:: FT
    b₁ #:: FT
    b₂ #:: FT
    b₃ #:: FT
    c₀ #:: FT
    c₁ #:: FT
    c₂ #:: FT
    d₀ #:: FT
    d₁ #:: FT
    k₀ #:: FT
    k₁ #:: FT
    k₂ #:: FT
    p₀ #:: FT
    v₀ #:: FT
    v₁ #:: FT
    v₂ #:: FT
    v₃ #:: FT
end

"""
    CarbonCoefficientParameters(a₀=0, a₁=0, a₂=0, a₃=0, a₄=0, a₅=0, b₀=0, b₁=0, b₂=0, b₃=0, c₀=0, c₁=0, c₂=0, d₀=0, d₁=0, k₀=0, k₁=0, k₂=0, p₀=0, v₀=0, v₁=0, v₂=0, v₃=0)

Create a `CarbonCoefficientParameters` object with the specified coefficients, converting them to the element type of the provided `grid`.

# Arguments
- `grid`: The grid object whose element type will be used for the coefficients.
- `a₀`, `a₁`, `a₂`, `a₃`, `a₄`, `a₅`: Coefficients for the `a` parameters (default is 0).
- `b₀`, `b₁`, `b₂`, `b₃`: Coefficients for the `b` parameters (default is 0).
- `c₀`, `c₁`, `c₂`: Coefficients for the `c` parameters (default is 0).
- `d₀`, `d₁`: Coefficients for the `d` parameters (default is 0).
- `k₀`, `k₁`, `k₂`: Coefficients for the `k` parameters (default is 0).
- `p₀`: Coefficient for the `p` parameter (default is 0).
- `v₀`, `v₁`, `v₂`, `v₃`: Coefficients for the `v` parameters (default is 0).

# Returns
- A `CarbonCoefficientParameters` object with the specified coefficients converted to the element type of the provided `grid`.
"""
function CarbonCoefficientParameters(;
    a₀ = 0.,
    a₁ = 0.,
    a₂ = 0.,
    a₃ = 0.,
    a₄ = 0.,
    a₅ = 0.,
    b₀ = 0.,
    b₁ = 0.,
    b₂ = 0.,
    b₃ = 0.,
    c₀ = 0.,
    c₁ = 0.,
    c₂ = 0.,
    d₀ = 0.,
    d₁ = 0.,
    k₀ = 0.,
    k₁ = 0.,
    k₂ = 0.,
    p₀ = 0.,
    v₀ = 0.,
    v₁ = 0.,
    v₂ = 0.,
    v₃ = 0.,
)
    return CarbonCoefficientParameters(
        a₀,
        a₁,
        a₂,
        a₃,
        a₄,
        a₅,
        b₀,
        b₁,
        b₂,
        b₃,
        c₀,
        c₁,
        c₂,
        d₀,
        d₁,
        k₀,
        k₁,
        k₂,
        p₀,
        v₀,
        v₁,
        v₂,
        v₃,
    )
end

Base.@kwdef struct CarbonSystemParameters #{CSP<:CarbonSolverParameters, CCP<:CarbonCoefficientParameters}
    Sᵒᵖᵗˢ        #:: CSP
    Pᴴ²⁰ˢʷ       #:: CCP
    Pᵘˢ          #:: CCP
    Pᴮᵀᴼᵀ        #:: CCP
    Pᶜᵃᵀᴼᵀ       #:: CCP
    Pᶠᵀᴼᵀ        #:: CCP
    Pˢᴼ⁴ᵀᴼᵀ      #:: CCP
    Pᵈⁱᶜₖₛₒₗₐ    # :: CCP
    Pᵈⁱᶜₖₚᵣₑ     #:: CCP
    Pᵈⁱᶜₖ₀       #:: CCP
    Pᵈⁱᶜₖ₁ᵣ₉₃    #:: CCP
    Pᵈⁱᶜₖ₂ᵣ₉₃    #:: CCP
    Pᵈⁱᶜₖ₁ₘ₉₅    #:: CCP
    Pᵈⁱᶜₖ₂ₘ₉₅    #:: CCP
    Pᵈⁱᶜₖ₁ₗ₀₀    # :: CCP
    Pᵈⁱᶜₖ₂ₗ₀₀    # :: CCP
    Pᴮₖ₁         #:: CCP
    Pᴴ²ᴼₖ₁       #:: CCP
    Pᴾᴼ⁴ₖ₁       #:: CCP
    Pᴾᴼ⁴ₖ₂       #:: CCP
    Pᴾᴼ⁴ₖ₃       #:: CCP
    Pˢⁱᵗₖ₁       #:: CCP
    Pᴴ²ˢₖ₁       #:: CCP
    Pᴺᴴ⁴ₖ₁       #:: CCP
    Pᴴᶠᵦ₁        #:: CCP
    Pᴴᶠₖ₁        #:: CCP
    Pᴴˢᴼ⁴ₖ₁      #:: CCP
    Pᶜᵃˡᶜⁱᵗᵉₛₚ   #:: CCP
    Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ #:: CCP
end
"""
    CarbonSystemParameters(kwargs...)

This function initializes and returns a `CarbonSystemParameters` object with various carbon coefficient parameters. 
The parameters are used in carbon system calculations and are initialized with default values, which can be overridden 
    by passing keyword arguments.

# Arguments
- `grid`: The grid on which the carbon system parameters are defined.
- `kwargs...`: Keyword arguments for various carbon coefficient parameters.

# Keyword Arguments
- `Sᵒᵖᵗˢ::CarbonSolverParameters`: Default is `CarbonSolverParameters(grid)`.
- `Pᴴ²⁰ˢʷ::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=1.0, a₁=-0.001005)`.
- `Pᵘˢ::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=0.019924)`.
- `Pᴮᵀᴼᵀ::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=0.000416, a₁=35.0, a₂=1.0)`.
- `Pᶜᵃᵀᴼᵀ::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=0.02127, a₁=40.078, a₂=1.80655)`.
- `Pᶠᵀᴼᵀ::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=6.8e-5, a₁=35.0)`.
- `Pˢᴼ⁴ᵀᴼᵀ::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=0.1400, a₁=96.062, a₂=1.80655)`.
- `Pᵈⁱᶜₖₛₒₗₐ::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=-162.8301, a₁=218.2968, a₂=90.9241, a₃=-1.47696, b₀=0.025695, b₁=-0.025225, b₂=0.0049867)`.
- `Pᵈⁱᶜₖₚᵣₑ::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=-1636.75, a₁=-12.0408, a₂=-0.0327957, a₃=3.16528e-5, b₀=57.7, b₁=-0.118, p₀=1.01325)`.
- `Pᵈⁱᶜₖ₀::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=-60.2409, a₁=93.4517, a₂=23.3585, b₀=0.023517, b₁=-0.023656, b₂=0.0047036)`.
- `Pᵈⁱᶜₖ₁ᵣ₉₃::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=2.83655, a₁=-2307.1266, a₂=-1.5529413, b₀=-0.20760841, b₁=-4.0484, b₂=0.08468345, b₃=-0.00654208, v₀=-25.5, v₁=-0.151, v₂=0.1271, k₀=-3.08e-3, k₁=-0.578e-3, k₂=0.0877e-3)`.
- `Pᵈⁱᶜₖ₂ᵣ₉₃::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=-9.226508, a₁=-3351.6106, a₂=-0.2005743, b₀=-0.106901773, b₁=-23.9722, b₂=0.1130822, b₃=-0.00846934, v₀=-15.82, v₁=0.321, v₂=-0.0219, k₀=1.13e-3, k₁=-0.314e-3, k₂=-0.1475e-3)`.
- `Pᵈⁱᶜₖ₁ₘ₉₅::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=2.18867, a₁=-2275.0360, a₂=-1.468591, b₀=-0.138681, b₁=-9.33291, b₂=0.0726483, b₃=-0.00574938, v₀=-25.5, v₁=-0.151, v₂=0.1271, k₀=-3.08e-3, k₁=-0.578e-3, k₂=0.0877e-3)`.
- `Pᵈⁱᶜₖ₂ₘ₉₅::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=-0.84226, a₁=-3741.1288, a₂=-1.437139, b₀=-0.128417, b₁=-24.41239, b₂=0.1195308, b₃=-0.00912840, v₀=-15.82, v₁=0.321, v₂=-0.0219, k₀=1.13e-3, k₁=-0.314e-3, k₂=-0.1475e-3)`.
- `Pᵈⁱᶜₖ₁ₗ₀₀::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=61.2172, a₁=-3633.86, a₂=-9.67770, b₀=0.011555, b₁=-0.0001152, v₀=-25.5, v₁=-0.151, v₂=0.1271, k₀=-3.08e-3, k₁=-0.578e-3, k₂=0.0877e-3)`.
- `Pᵈⁱᶜₖ₂ₗ₀₀::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=-25.9290, a₁=-471.78, a₂=3.16967, b₀=0.01781, b₁=-0.0001122, v₀=-15.82, v₁=0.321, v₂=-0.0219, k₀=1.13e-3, k₁=-0.314e-3, k₂=-0.1475e-3)`.
- `Pᴮₖ₁::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=-8966.90, a₁=-2890.53, a₂=-77.942, a₃=1.728, a₄=-0.0996, b₀=148.0248, b₁=137.1942, b₂=1.62142, c₀=-24.4344, c₁=-25.085, c₂=-0.2474, d₀=0.053105, v₀=-29.48, v₁=0.295, v₂=0.1622, v₃=-0.002608, k₀=-2.84e-3, k₁=0.354e-3)`.
- `Pᴴ²ᴼₖ₁::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=148.9802, a₁=-13847.26, a₂=-23.6521, b₀=-5.977, b₁=118.67, b₂=1.0495, c₀=-0.01615, v₀=-20.02, v₁=0.1119, v₂=-0.1409E-02, k₀=-5.13e-3, k₁=0.0794e-3)`.
- `Pᴾᴼ⁴ₖ₁::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=115.54, a₁=-4576.752, a₂=-18.453, b₀=0.69171, b₁=-106.736, b₂=-0.01844, b₃=-0.65643, v₀=-14.51, v₁=0.1211, v₂=-0.321E-03, k₀=-2.67e-3, k₁=0.0427e-3)`.
- `Pᴾᴼ⁴ₖ₂::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=172.1033, a₁=-8814.715, a₂=-27.927, b₀=1.3566, b₁=-160.340, b₂=-0.05778, b₃=0.37335, v₀=-23.12, v₁=0.1758, v₂=-0.002647, k₀=-5.15e-3, k₁=0.09e-3)`.
- `Pᴾᴼ⁴ₖ₃::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=-18.126, a₁=-3070.75, a₂=2.81197, a₃=17.27039, a₄=-0.09984, a₅=-44.99486, v₀=-26.57, v₁=0.2020, v₂=-3.042e-3, k₀=-4.08e-3, k₁=0.0714e-3)`.
- `Pˢⁱᵗₖ₁::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=117.40, a₁=-8904.2, a₂=-19.334, b₀=3.5913, b₁=-458.79, b₂=-1.5998, b₃=188.74, c₀=0.07871, c₁=-12.1652, v₀=-29.48, v₁=0.0, v₂=0.1622, v₃=-0.002608, k₀=-2.84e-3, k₁=0.354e-3)`.
- `Pᴴ²ˢₖ₁::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=225.838, a₁=-13275.3, a₂=-34.6435, a₃=0.3449, a₄=-0.0274, v₀=-14.80, v₁=0.0020, v₂=-0.400E-03, k₀=2.89e-3, k₁=0.054e-3)`.
- `Pᴺᴴ⁴ₖ₁::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=-0.25444, a₁=-6285.33, a₂=0.0001635, b₀=0.46532, b₁=-123.7184, b₂=-0.01992, b₃=3.17556, v₀=-26.43, v₁=0.0889, v₂=-0.905E-03, k₀=-5.03E-03, k₁=0.0814E-03)`.
- `Pᴴᶠᵦ₁::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=12.641, a₁=-1590.2, a₂=-1.525, v₀=-9.78, v₁=-0.0090, v₂=-0.942E-03, k₀=-3.91e-3, k₁=0.054e-3)`.
- `Pᴴᶠₖ₁::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=-9.68, a₁=874.0, a₂=0.111, v₀=-9.78, v₁=-0.0090, v₂=-0.942E-3, k₀=-3.91e-3, k₁=0.054e-3)`.
- `Pᴴˢᴼ⁴ₖ₁::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=141.328, a₁=-4276.1, a₂=-23.093, b₀=324.57, b₁=-13856., b₂=-47.986, c₀=-771.54, c₁=35474., c₂=114.723, d₀=-2698., d₁=1776., v₀=-18.03, v₁=0.0466, v₂=0.316E-03, k₀=-4.53e-3, k₁=0.0900e-3)`.
- `Pᶜᵃˡᶜⁱᵗᵉₛₚ::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=-171.9065, a₁=-0.077993, a₂=2839.319, a₃=71.595, b₀=-0.77712, b₁=0.0028426, b₂=178.34, c₀=-0.07711, d₀=0.0041249, v₀=-48.76, v₁=0.5304, k₀=-11.76e-3, k₁=0.3692e-3)`.
- `Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ::CarbonCoefficientParameters`: Default is `CarbonCoefficientParameters(a₀=-171.945, a₁=-0.077993, a₂=2903.293, a₃=71.595, b₀=-0.068393, b₁=0.0017276, b₂=88.135, c₀=-0.10018, d₀=0.0059415, v₀=-48.76, v₁=0.5304, v₂=2.8, k₀=-11.76e-3, k₁=0.3692e-3)`.

# Returns
- `CarbonSystemParameters`: An object containing all the initialized carbon coefficient parameters.
"""
function CarbonSystemParameters(;
        Sᵒᵖᵗˢ :: CarbonSolverParameters = CarbonSolverParameters(),
        Pᴴ²⁰ˢʷ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
            a₀ =   1.0,
            a₁ = - 0.001005,
        ),
        Pᵘˢ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
            a₀ =   0.019924,
        ),
        Pᴮᵀᴼᵀ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
            a₀ =   0.000416,
            a₁ =   35.0,
            a₂ =   1.0,
        ),
        Pᶜᵃᵀᴼᵀ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
            a₀ =   0.02127,
            a₁ =   40.078,
            a₂ =   1.80655,
        ),
        Pᶠᵀᴼᵀ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
            a₀ =   6.8e-5,
            a₁ =   35.0,
        ),
        Pˢᴼ⁴ᵀᴼᵀ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
            a₀ =   0.1400,
            a₁ =   96.062,
            a₂ =   1.80655,
        ),
        Pᵈⁱᶜₖₛₒₗₐ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
            a₀ = -162.8301,
            a₁ =  218.2968,
            a₂ =   90.9241,
            a₃ = -  1.47696,
            b₀ =   0.025695,
            b₁ = - 0.025225,
            b₂ =   0.0049867,
        ),
        Pᵈⁱᶜₖₚᵣₑ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
            a₀ = -1636.75,
            a₁ = -  12.0408,
            a₂ = -   0.0327957, 
            a₃ =     3.16528e-5,
            b₀ =    57.7,
            b₁ = -   0.118,
            p₀ =     1.01325, # p_bar_oneatmosphere, Handbook (2007)
        ),
        Pᵈⁱᶜₖ₀ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
            a₀ = -60.2409,
            a₁ =  93.4517,
            a₂ =  23.3585,
            b₀ =   0.023517,
            b₁ = - 0.023656,
            b₂ =   0.0047036,
        ),
        Pᵈⁱᶜₖ₁ᵣ₉₃ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
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
        ),
        Pᵈⁱᶜₖ₂ᵣ₉₃ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
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
        ),
        Pᵈⁱᶜₖ₁ₘ₉₅ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
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
        ),
        Pᵈⁱᶜₖ₂ₘ₉₅ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
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
        ),
        Pᵈⁱᶜₖ₁ₗ₀₀ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
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
        ),
        Pᵈⁱᶜₖ₂ₗ₀₀ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
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
        ),
        Pᴮₖ₁ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
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
        ),
        Pᴴ²ᴼₖ₁ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
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
            k₁ =      0.0794e-3,
        ),
        Pᴾᴼ⁴ₖ₁ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
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
            k₁ =     0.0427e-3,
        ),
        Pᴾᴼ⁴ₖ₂ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
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
            k₁ =      0.09e-3,
        ),
        Pᴾᴼ⁴ₖ₃ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
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
            k₁ =      0.0714e-3,
        ),
        Pˢⁱᵗₖ₁ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
            a₀ =    117.40,
            a₁ = - 8904.2,
            a₂ = -   19.334,
            b₀ =      3.5913,
            b₁ = -  458.79,
            b₂ = -    1.5998,
            b₃ =    188.74,
            c₀ =      0.07871,
            c₁ = -   12.1652,
            v₀ = -   29.48,
            v₁ =      0.0,
            v₂ =      0.1622,
            v₃ = -    0.002608,
            k₀ = -    2.84e-3,
            k₁ =      0.354e-3,
        ),
        Pᴴ²ˢₖ₁ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
            a₀ =     225.838,
            a₁ = - 13275.3,
            a₂ = -    34.6435,
            a₃ =       0.3449,
            a₄ = -     0.0274,
            v₀ = -    14.80,
            v₁ =       0.0020,
            v₂ =  -    0.400E-03,
            k₀ =       2.89e-3,
            k₁ =       0.054e-3,
        ),
        Pᴺᴴ⁴ₖ₁ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
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
            k₁ =      0.0814E-03,
        ),
        Pᴴᶠᵦ₁ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
            a₀ =     12.641,
            a₁ = - 1590.2,
            a₂ = -    1.525,
            v₀ = -    9.78,
            v₁ = -    0.0090,
            v₂ = -    0.942E-03,
            k₀ = -    3.91e-3,
            k₁ =      0.054e-3,
        ),
        Pᴴᶠₖ₁ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
            a₀ = -  9.68,
            a₁ =  874.0,
            a₂ =    0.111,
            v₀ = -  9.78,
            v₁ = -  0.0090,
            v₂ = -  0.942E-3,
            k₀ = -  3.91e-3,
            k₁ =    0.054e-3,
        ),
        Pᴴˢᴼ⁴ₖ₁ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
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
            k₁ =       0.0900e-3,
        ),
        Pᶜᵃˡᶜⁱᵗᵉₛₚ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
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
            k₁ =     0.3692e-3,
        ),
        Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ :: CarbonCoefficientParameters = CarbonCoefficientParameters( # 
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
            k₁ =     0.3692e-3,
        )
)
    return CarbonSystemParameters(
        Sᵒᵖᵗˢ,
        Pᴴ²⁰ˢʷ, 
        Pᵘˢ,
        Pᴮᵀᴼᵀ,
        Pᶜᵃᵀᴼᵀ,
        Pᶠᵀᴼᵀ,
        Pˢᴼ⁴ᵀᴼᵀ,
        Pᵈⁱᶜₖₛₒₗₐ,
        Pᵈⁱᶜₖₚᵣₑ,
        Pᵈⁱᶜₖ₀,
        Pᵈⁱᶜₖ₁ᵣ₉₃,
        Pᵈⁱᶜₖ₂ᵣ₉₃,
        Pᵈⁱᶜₖ₁ₘ₉₅,
        Pᵈⁱᶜₖ₂ₘ₉₅,
        Pᵈⁱᶜₖ₁ₗ₀₀,
        Pᵈⁱᶜₖ₂ₗ₀₀,
        Pᴮₖ₁,
        Pᴴ²ᴼₖ₁,
        Pᴾᴼ⁴ₖ₁,
        Pᴾᴼ⁴ₖ₂,
        Pᴾᴼ⁴ₖ₃,
        Pˢⁱᵗₖ₁,
        Pᴴ²ˢₖ₁,
        Pᴺᴴ⁴ₖ₁,
        Pᴴᶠᵦ₁,
        Pᴴᶠₖ₁,
        Pᴴˢᴼ⁴ₖ₁,
        Pᶜᵃˡᶜⁱᵗᵉₛₚ,
        Pᵃʳᵃᵍᵒⁿⁱᵗᵉₛₚ,
    )
end