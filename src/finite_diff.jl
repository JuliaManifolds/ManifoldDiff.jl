
"""
    FiniteDiffBackend <: AbstractDiffBackend

A type to specify / use differentiation backend based on FiniteDiff.jl package.

# Constructor

    FiniteDiffBackend(method::Val{Symbol} = Val{:central})
"""
struct FiniteDiffBackend{TM<:Val} <: AbstractDiffBackend
    method::TM
end

FiniteDiffBackend() = FiniteDiffBackend(Val(:central))
