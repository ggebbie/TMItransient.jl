# TMItransient

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ggebbie.github.io/TMItransient.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ggebbie.github.io/TMItransient.jl/dev)
[![Build Status](https://github.com/ggebbie/TMItransient.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ggebbie/TMItransient.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ggebbie/TMItransient.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ggebbie/TMItransient.jl)

Transient simulations using the TMI circulation matrix. This Julia package builds upon 
the package [TMI.jl](https://github.com/ggebbie/TMI.jl).

# Package setup

This package was setup with PkgTemplates.jl. The starting code was

```
using PkgTemplates

t = Template(; user="ggebbie", dir="~/projects", authors="G Jake Gebbie", julia=v"1.7", plugins=[ License(; name="MIT"), Git(; manifest=true, ssh=true), GitHubActions(; x86=false, extra_versions=["1.7","nightly"]), Codecov(), Documenter{GitHubActions}(), Develop(), ], )

t("TMItransient.jl")
```

Follow [notes](https://m3g.github.io/JuliaNotes.jl/stable/publish_docs/) to make documenter key and deploy key 
`import DocumenterTools`\
`DocumenterTools.genkeys()`\
`DocumenterTools.genkeys(user="ggebbie", repo="TMItransient.jl")` # optionally (I think)\
Note: must call second key "DOCUMENTER_KEY"

- optional: add argument to deploydocs in docs/make.jl "devbranch="main" or "numerics" etc. Make a /dev version of docs. Will it make a stable version when a release is made?

- to build docs manually, try julia --project=docs docs/make.jl when I did it locally, I activated TMI project, then include("make.jl") and it worked locally

- previously I made a gh-pages branch following, but it happened automatically this time [[https://coderwall.com/p/0n3soa/create-a-disconnected-git-branch][instructions for creating a disconnected git branch]]

- You can also update the docs just by uploading a new tag, with:

`git tag -a v0.1.0+doc2 -m "v0.1.0"`\
`git push --tag`
