module ManifoldDiffForwardDiffExt

using ManifoldDiff
using ManifoldDiff: ForwardDiffBackend
using ForwardDiff

function ManifoldDiff._derivative(f, p, ::ForwardDiffBackend)
    return ForwardDiff.derivative(f, p)
end

function ManifoldDiff._derivative!(f, X, t, ::ForwardDiffBackend)
    return ForwardDiff.derivative!(X, f, t)
end

function ManifoldDiff._gradient(f, p, ::ForwardDiffBackend)
    return ForwardDiff.gradient(f, p)
end

function ManifoldDiff._gradient!(f, X, t, ::ForwardDiffBackend)
    return ForwardDiff.gradient!(X, f, t)
end

function ManifoldDiff._jacobian(f, p, ::ForwardDiffBackend)
    return ForwardDiff.jacobian(f, p)
end

function ManifoldDiff._hessian(f, p, ::ForwardDiffBackend)
    return ForwardDiff.hessian(f, p)
end

end
