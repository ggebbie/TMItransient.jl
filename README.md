# TMItransient

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ggebbie.github.io/TMItransient.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ggebbie.github.io/TMItransient.jl/dev)
[![Build Status](https://github.com/ggebbie/TMItransient.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ggebbie/TMItransient.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ggebbie/TMItransient.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ggebbie/TMItransient.jl)

# Package setup

This package was setup with PkgTemplates.jl. The starting code was

`using PkgTemplates

t = Template(; user="ggebbie", dir="~/projects", authors="G Jake Gebbie", julia=v"1.7", plugins=[ License(; name="MIT"), Git(; manifest=true, ssh=true), GitHubActions(; x86=false, extra_versions=["1.7","nightly"]), Codecov(), Documenter{GitHubActions}(), Develop(), ], )

t("TMItransient.jl")`
