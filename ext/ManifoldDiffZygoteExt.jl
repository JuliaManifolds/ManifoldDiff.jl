module ManifoldDiffZygoteExt

using ManifoldDiff
using ManifoldDiff: ZygoteDiffBackend
using Zygote

function ManifoldDiff._gradient(f, p, ::ZygoteDiffBackend)
    return Zygote.gradient(f, p)[1]
end

function ManifoldDiff._gradient!(f, X, p, ::ZygoteDiffBackend)
    return copyto!(X, Zygote.gradient(f, p)[1])
end


end
