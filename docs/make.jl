using Plots, Manifolds, ManifoldsBase, ManifoldDiff, Documenter, PyPlot
ENV["GKSwstype"] = "100"

makedocs(
    # for development, we disable prettyurls
    format = Documenter.HTML(prettyurls = false, assets = ["assets/favicon.ico"]),
    modules = [ManifoldDiffDev],
    authors = "Seth Axen, Mateusz Baran, Ronny Bergmann, and contributors.",
    sitename = "ManifoldDiff.jl",
    pages = ["Home" => "index.md", "Library" => "library.md"],
)
deploydocs(
    repo = "github.com/JuliaManifolds/ManifoldDiff.jl.git",
    push_preview = true,
    devbranch = "main",
)
