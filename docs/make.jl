using Documenter
using QuantumInputOutput, SecondQuantizedAlgebra

ENV["GKSwstype"] = "100" # enable headless mode for GR to suppress warnings when plotting

include("make_md_examples.jl")

pages = [
    "index.md",
    "theory.md",
    "tutorial.md",
    "implementation.md",
    "api.md",
    "Examples" => [
        "examples/01-1_cavity-scattering__PRL2019_123-123604_fig2-fig3.md",
        "examples/01-2_stimulated-emission__PRL2019_123-123604_fig4.md",
        "examples/02-1_cavity-phase-noise__PRA2020_102- 023717_fig2.md"
    ],
]

using Pkg
status = sprint(io -> Pkg.status("SecondQuantizedAlgebra"; io = io))
# version = match(r"(v[0-9].[0-9]+.[0-9]+)", status)[1]
gh_moi = Documenter.Remotes.GitHub("qojulia", "SecondQuantizedAlgebra.jl")
# remotes = Dict(pkgdir(SecondQuantizedAlgebra) => (gh_moi, version))

makedocs(
    sitename = "QuantumInputOutput.jl",
    modules = [QuantumInputOutput],#, SecondQuantizedAlgebra],
    pages = pages,
    # remotes = remotes,
    checkdocs = :exports,
    format = Documenter.HTML(
        mathengine = MathJax(),
        footer = "[**Back to GitHub**](https://github.com/ChristophHotter/QuantumInputOutput.jl)",
        example_size_threshold = 800 * 2^10,
        size_threshold_warn = 400 * 2^10,
        size_threshold = 600 * 2^10,
    ),
)

deploydocs(repo = "github.com/ChristophHotter/QuantumInputOutput.jl", push_preview = true)
