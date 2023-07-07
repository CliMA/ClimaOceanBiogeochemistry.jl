module ClimaOceanBiogeochemistry

# solvers for ocean pH and pCO2
include("CarbonSystemSolvers.jl")

# diffferent kinds of BGC model
include("nutrients_plankton_bacteria_detritus.jl")
include("carbon_alkalinity_nutrients.jl")

end # module ClimaOceanBiogeochemistry
