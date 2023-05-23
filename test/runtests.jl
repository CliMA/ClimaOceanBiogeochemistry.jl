using Test
using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients, NutrientsPlanktonBacteriaDetritus

@testset "CarbonSolvers" begin
    include("test_carbon_solvers.jl")
end

@test NutrientsPlanktonBacteriaDetritus() isa NutrientsPlanktonBacteriaDetritus
@test CarbonAlkalinityNutrients() isa CarbonAlkalinityNutrients
