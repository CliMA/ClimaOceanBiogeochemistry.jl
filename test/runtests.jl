using Test

@test 1 == 1

@testset "CarbonSolvers" begin
    include("test_carbon_solvers.jl")
end
