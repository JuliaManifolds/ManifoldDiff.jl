#!/usr/bin/env julia
#
#

if "--help" ∈ ARGS
    println(
        """
     docs/make.jl

Render the `ManifoldDiff.jl` documentation with optional arguments

Arguments
* `--help`         - print this help and exit without rendering the documentation
* `--prettyurls`   – toggle the prettyurls part to true (which is otherwise only true on CI)
""",
    )
    exit(0)
end
# (a) if docs is not the current active environment, switch to it
# (from https://github.com/JuliaIO/HDF5.jl/pull/1020/) 
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    using Pkg
    Pkg.activate(@__DIR__)
    Pkg.develop(PackageSpec(; path = (@__DIR__) * "/../"))
    Pkg.resolve()
    Pkg.instantiate()
end

using ManifoldsBase, ManifoldDiff
using Documenter, DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "references.bib"); style = :alpha)
makedocs(;
    format = Documenter.HTML(;
        prettyurls = (get(ENV, "CI", nothing) == "true") || ("--prettyurls" ∈ ARGS),
        assets = ["assets/favicon.ico", "assets/citations.css"],
    ),
    modules = [ManifoldDiff],
    authors = "Seth Axen, Mateusz Baran, Ronny Bergmann, and contributors.",
    sitename = "ManifoldDiff.jl",
    pages = [
        "Home" => "index.md",
        "Usage" => "basic_usage.md",
        "Backends" => "backends.md",
        "Library of functions" => "library.md",
        "Internals" => "internals.md",
        "Literature" => "references.md",
    ],
    plugins = [bib],
    warnonly = [:autodocs_block],
)
deploydocs(
    repo = "github.com/JuliaManifolds/ManifoldDiff.jl.git",
    push_preview = true,
    devbranch = "main",
)
