module ManifoldDiffFiniteDiffExt

using ManifoldDiff
using ManifoldDiff: FiniteDiffBackend
using FiniteDiff

function ManifoldDiff._derivative(f, p, ::FiniteDiffBackend{Method}) where {Method}
    return FiniteDiff.finite_difference_derivative(f, p, Method)
end

function ManifoldDiff._gradient(f, p, ::FiniteDiffBackend{Method}) where {Method}
    return FiniteDiff.finite_difference_gradient(f, p, Method)
end

function ManifoldDiff._gradient!(f, X, p, ::FiniteDiffBackend{Method}) where {Method}
    return FiniteDiff.finite_difference_gradient!(X, f, p, Method)
end

function ManifoldDiff._jacobian(f, p, ::FiniteDiffBackend{Method}) where {Method}
    return FiniteDiff.finite_difference_jacobian(f, p, Method)
end

function ManifoldDiff._jacobian!(f, X, p, ::FiniteDiffBackend{Method}) where {Method}
    return FiniteDiff.finite_difference_jacobian!(X, f, p, Method)
end


end
