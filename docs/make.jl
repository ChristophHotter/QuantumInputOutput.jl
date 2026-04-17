using Documenter
using QuantumInputOutput, SecondQuantizedAlgebra

using Plots
default(; fmt = :png)

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
        "examples/02-1_cavity-phase-noise__PRA2020_102- 023717_fig2.md",
        "examples/02-3_mode-entanglement__PRA2020_102-023717_fig4.md",
        "examples/03-1_beam-combiner__PRA2023_107-023715_fig2-fig3.md",
        "examples/04-1_two-sided-cavity_with-atom_coh-drive.md",
        "examples/04-2_two-sided-cavity_with-atom_coh-drive__cumulants.md",
        "examples/05-1_N-QDs_bidirectional-waveguide_coherent-pulse.md",
        "examples/05-2_N-QDs_bidirectional-waveguide_quantum-pulse_qo.md",
        "examples/05-3_N-QDs_bidirectional-waveguide_feedback-reduction.md",
        "examples/06-1_interaction-picture__PRA2023_107-013706_fig2.md",
        "examples/07-1_beamsplitter_loss__quantum-pulse.md",
        "examples/07-2_hong-ou-mandel__quantum-pulse.md",
        "examples/08-1_pulse-delay__simple.md",
        "examples/09-1_coherent-feedback-squeezing__Gough-Wildfeuer-2009.md",
    ],
]

makedocs(
    sitename = "QuantumInputOutput.jl",
    modules = QuantumInputOutput,
    format = Documenter.HTML(;
        canonical = "https://qojulia.github.io/QuantumInputOutput.jl",
    ),
    pages = pages,
    clean = true,
    linkcheck = false,
    warnonly = :missing_docs,
    draft = false,#,(!CI),
    doctest = false,  # We test it in the CI, no need to run it here
    checkdocs = :exports,
)

deploydocs(
    repo = "github.com/qojulia/QuantumInputOutput.jl",
    devbranch = "main",
    target = "build",
    branch = "gh-pages",
    push_preview = true,
)
