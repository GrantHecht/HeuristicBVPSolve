using HeuristicBVPSolve
using Documenter

DocMeta.setdocmeta!(HeuristicBVPSolve, :DocTestSetup, :(using HeuristicBVPSolve); recursive=true)

makedocs(;
    modules=[HeuristicBVPSolve],
    authors="Grant Hecht",
    repo="https://github.com/GrantHecht/HeuristicBVPSolve.jl/blob/{commit}{path}#{line}",
    sitename="HeuristicBVPSolve.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://GrantHecht.github.io/HeuristicBVPSolve.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/GrantHecht/HeuristicBVPSolve.jl",
    devbranch="main",
)
