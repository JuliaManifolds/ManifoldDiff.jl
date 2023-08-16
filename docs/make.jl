using Manifolds, ManifoldsBase, ManifoldDiff, Documenter
using ForwardDiff, ReverseDiff, FiniteDifferences, Zygote

makedocs(
    # for development, we disable prettyurls
    format = Documenter.HTML(prettyurls = false, assets = ["assets/favicon.ico"]),
    modules = [ManifoldDiff],
    authors = "Seth Axen, Mateusz Baran, Ronny Bergmann, and contributors.",
    sitename = "ManifoldDiff.jl",
    pages = ["Home" => "index.md", "Library" => "library.md"],
)
deploydocs(
    repo = "github.com/JuliaManifolds/ManifoldDiff.jl.git",
    push_preview = true,
    devbranch = "main",
)
