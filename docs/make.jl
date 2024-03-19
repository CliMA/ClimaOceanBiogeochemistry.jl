using
  Documenter,
  Literate,
  ClimaOceanBiogeochemistry,
  ClimaOceanBiogeochemistry.CarbonSystemSolvers

#####
##### Generate examples
#####

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/literated")

to_be_literated = [
    "single_column_nutrients_plankton_bacteria_detritus.jl",
]
  
for file in to_be_literated
    filepath = joinpath(EXAMPLES_DIR, file)
    Literate.markdown(filepath, OUTPUT_DIR;
                      flavor = Literate.DocumenterFlavor())
end

#####
##### Build and deploy docs
#####

MiB = 2^20

format = Documenter.HTML(collapselevel = 2,
                         prettyurls = get(ENV, "CI", nothing) == "true",
                         size_threshold = 1MiB,
                         canonical = "https://clima.github.io/ClimaOceanBiogeochemistry.jl/dev/")

pages = [
    "Home" => "index.md",
    "Examples" => [
        "Single column nutrients, plankton, bacteria, detritus" => "literated/single_column_nutrients_plankton_bacteria_detritus.md",
    ],
    "Library" => [ 
        "Contents" => "library/outline.md",
        "Public" => "library/public.md",
        "Private" => "library/internals.md",
        "Function index" => "library/function_index.md",
    ],
]

makedocs(sitename = "ClimaOceanBiogeochemistry.jl",
         modules = [ClimaOceanBiogeochemistry, 
                    ClimaOceanBiogeochemistry.CarbonSystemSolvers,
                    ClimaOceanBiogeochemistry.CarbonSystemSolvers.UniversalRobustCarbonSolver,
                    ClimaOceanBiogeochemistry.CarbonSystemSolvers.AlkalinityCorrectionCarbonSolver,
                    ClimaOceanBiogeochemistry.CarbonSystemSolvers.DirectCubicCarbonSolver],
         format = format,
         pages = pages,
         doctest = true,
         warnonly = [:cross_references],
         clean = true,
         checkdocs = :exports)


@info "Clean up temporary .jld2/.nc files created by doctests..."

"""
    recursive_find(directory, pattern)

Return list of filepaths within `directory` that contains the `pattern::Regex`.
"""
recursive_find(directory, pattern) =
    mapreduce(vcat, walkdir(directory)) do (root, dirs, files)
        joinpath.(root, filter(contains(pattern), files))
    end

files = []
for pattern in [r"\.jld2", r"\.nc"]
    global files = vcat(files, recursive_find(@__DIR__, pattern))
end

for file in files
    rm(file)
end


withenv("GITHUB_REPOSITORY" => "CliMA/ClimaOceanBiogeochemistry.jl") do
    deploydocs(repo = "github.com/CliMA/ClimaOceanBiogeochemistry.jl.git",
               versions = ["stable" => "v^", "dev" => "dev", "v#.#.#"],
               forcepush = true,
               devbranch = "main",
               push_preview = true)
end
