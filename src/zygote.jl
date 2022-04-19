struct ZygoteDiffBackend <: AbstractDiffBackend end

function ManifoldDiff._gradient(f, p, ::ZygoteDiffBackend)
    return Zygote.gradient(f, p)[1]
end

function ManifoldDiff._gradient!(f, X, p, ::ZygoteDiffBackend)
    return copyto!(X, Zygote.gradient(f, p)[1])
end

if default_differential_backend() === NoneDiffBackend()
    set_default_differential_backend!(ZygoteDiffBackend())
end
