pushfirst!(LOAD_PATH, joinpath(@__DIR__, "..")) # add ClimaOceanBiogeochemistry to environment stack

using
  Documenter,
  Glob,
  Literate,
  ClimaOceanBiogeochemistry

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
    Literate.markdown(filepath, OUTPUT_DIR; flavor = Literate.DocumenterFlavor())
end

#####
##### Build and deploy docs
#####

# Set up a timer to print a space ' ' every 240 seconds. This is to avoid CI
# timing out when building demanding Literate.jl examples.
Timer(t -> println(" "), 0, interval=240)

format = Documenter.HTML(collapselevel = 2,
                         prettyurls = get(ENV, "CI", nothing) == "true",
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
         modules = [ClimaOceanBiogeochemistry],
         format = format,
         pages = pages,
         doctest = true,
         strict = true,
         clean = true,
         checkdocs = :exports)

withenv("GITHUB_REPOSITORY" => "CliMA/ClimaOceanBiogeochemistry.jl") do
    deploydocs(repo = "github.com/CliMA/ClimaOceanBiogeochemistry.jl.git",
               versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
               forcepush = true,
               devbranch = "main",
               push_preview = true)
end

