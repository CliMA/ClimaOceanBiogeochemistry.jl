using Test
using Oceananigans
using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients, NutrientsPlanktonBacteriaDetritus

@testset "CarbonSolvers" begin
    include("test_carbon_solvers.jl")
end

@testset "NutrientsPlanktonBacteriaDetritus" begin
    @test NutrientsPlanktonBacteriaDetritus() isa NutrientsPlanktonBacteriaDetritus
end

@testset "CarbonAlkalinityNutrients" begin
    @test CarbonAlkalinityNutrients() isa CarbonAlkalinityNutrients

    grid = RectilinearGrid(size = 64, z = (-256, 0), topology = (Flat, Flat, Bounded))
    model = HydrostaticFreeSurfaceModel(; grid, biogeochemistry = CarbonAlkalinityNutrients())

    @test :PO₄ ∈ keys(model.tracers)

    simulation = Simulation(model, Δt=10minutes, stop_iteration=3)

    @test try 
        run!(simulation)
        true
    catch
        false
    end
end

