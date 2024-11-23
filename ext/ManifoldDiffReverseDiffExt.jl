module ManifoldDiffReverseDiffExt

using ManifoldDiff
using ManifoldDiff: ReverseDiffBackend
using ReverseDiff

function ManifoldDiff._gradient(f, p, ::ReverseDiffBackend)
    return ReverseDiff.gradient(f, p)
end

function ManifoldDiff._gradient!(f, X, p, ::ReverseDiffBackend)
    return ReverseDiff.gradient!(X, f, p)
end

end
