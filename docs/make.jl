using TMItransient
using Documenter

DocMeta.setdocmeta!(TMItransient, :DocTestSetup, :(using TMItransient); recursive=true)

makedocs(;
    modules=[TMItransient],
    authors="G Jake Gebbie",
    repo="https://github.com/ggebbie/TMItransient.jl/blob/{commit}{path}#{line}",
    sitename="TMItransient.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ggebbie.github.io/TMItransient.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ggebbie/TMItransient.jl",
    devbranch="main",
)
