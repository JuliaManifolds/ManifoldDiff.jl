module ManifoldDiffReverseDiffExt

if isdefined(Base, :get_extension)
    using ManifoldDiff
    using ManifoldDiff: ReverseDiffBackend
    using ReverseDiff
else
    # imports need to be relative for Requires.jl-based workflows:
    # https://github.com/JuliaArrays/ArrayInterface.jl/pull/387
    using ..ManifoldDiff
    using ..ManifoldDiff: ReverseDiffBackend
    using ..ReverseDiff
end

function ManifoldDiff._gradient(f, p, ::ReverseDiffBackend)
    return ReverseDiff.gradient(f, p)
end

function ManifoldDiff._gradient!(f, X, p, ::ReverseDiffBackend)
    return ReverseDiff.gradient!(X, f, p)
end

end
