"""
    FiniteDifferencesBackend(method::FiniteDifferenceMethod = central_fdm(5, 1))

Differentiation backend based on the FiniteDifferences.jl package.
"""
struct FiniteDifferencesBackend{TM} <: AbstractDiffBackend
    method::TM
end
