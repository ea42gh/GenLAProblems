using GenLAProblems
using Documenter

DocMeta.setdocmeta!(GenLAProblems, :DocTestSetup, :(using GenLAProblems); recursive=true)

makedocs(;
    modules=[GenLAProblems],
    authors="ea42_github@mail.com",
    sitename="GenLAProblems.jl",
    format=Documenter.HTML(;
        canonical="https://ea42gh.github.io/GenLAProblems.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ea42gh/GenLAProblems.jl",
    devbranch="main",
)
