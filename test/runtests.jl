using Test
using Oceananigans
using Oceananigans.Units
using ClimaOceanBiogeochemistry: CarbonAlkalinityNutrients, NutrientsPlanktonBacteriaDetritus

@testset "CarbonSolvers" begin
    include("test_carbon_solvers.jl")
end

@testset "NutrientsPlanktonBacteriaDetritus" begin
    grid = RectilinearGrid(size=(4, 6, 8), extent=(1, 2, 3))
    @test NutrientsPlanktonBacteriaDetritus(grid) isa NutrientsPlanktonBacteriaDetritus
end

@testset "CarbonAlkalinityNutrients" begin
    grid = RectilinearGrid(size = 64, z = (-256, 0), topology = (Flat, Flat, Bounded))
    @test CarbonAlkalinityNutrients(; grid) isa CarbonAlkalinityNutrients

    model = HydrostaticFreeSurfaceModel(; grid, biogeochemistry = CarbonAlkalinityNutrients(; grid))

    @test :PO₄ ∈ keys(model.tracers)

    simulation = Simulation(model, Δt=10minutes, stop_iteration=3)

    @test try 
        run!(simulation)
        true
    catch
        false
    end
end
