module ManifoldDiffZygoteExt

if isdefined(Base, :get_extension)
    using ManifoldDiff
    using ManifoldDiff: ZygoteDiffBackend
    using Zygote
else
    # imports need to be relative for Requires.jl-based workflows:
    # https://github.com/JuliaArrays/ArrayInterface.jl/pull/387
    using ..ManifoldDiff
    using ..ManifoldDiff: ZygoteDiffBackend
    using ..Zygote
end

function ManifoldDiff._gradient(f, p, ::ZygoteDiffBackend)
    return Zygote.gradient(f, p)[1]
end

function ManifoldDiff._gradient!(f, X, p, ::ZygoteDiffBackend)
    return copyto!(X, Zygote.gradient(f, p)[1])
end


end
