# (a) if docs is not the current active environment, switch to it
# (from https://github.com/JuliaIO/HDF5.jl/pull/1020/) 
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    using Pkg
    Pkg.activate(@__DIR__)
    Pkg.develop(PackageSpec(; path = (@__DIR__) * "/../"))
    Pkg.resolve()
    Pkg.instantiate()
end

using ManifoldDiff, ManifoldsBase, ManifoldDiff
using Documenter, DocumenterCitations
using FiniteDiff, ForwardDiff, ReverseDiff, FiniteDifferences, Zygote

bib = CitationBibliography(joinpath(@__DIR__, "src", "references.bib"); style = :alpha)
makedocs(;
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = ["assets/favicon.ico"],
    ),
    modules = [
        ManifoldDiff,
        isdefined(Base, :get_extension) ?
        Base.get_extension(ManifoldDiff, :ManifoldDiffFiniteDiffExt) :
        ManifoldDiff.ManifoldDiffFiniteDiffExt,
        isdefined(Base, :get_extension) ?
        Base.get_extension(ManifoldDiff, :ManifoldDiffFiniteDifferencesExt) :
        ManifoldDiff.ManifoldDiffFiniteDifferencesExt,
        isdefined(Base, :get_extension) ?
        Base.get_extension(ManifoldDiff, :ManifoldDiffForwardDiffExt) :
        ManifoldDiff.ManifoldDiffForwardDiffExt,
        isdefined(Base, :get_extension) ?
        Base.get_extension(ManifoldDiff, :ManifoldDiffReverseDiffExt) :
        ManifoldDiff.ManifoldDiffReverseDiffExt,
        isdefined(Base, :get_extension) ?
        Base.get_extension(ManifoldDiff, :ManifoldDiffZygoteExt) :
        ManifoldDiff.ManifoldDiffZygoteExt,
    ],
    authors = "Seth Axen, Mateusz Baran, Ronny Bergmann, and contributors.",
    sitename = "ManifoldDiff.jl",
    pages = [
        "Home" => "index.md",
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
