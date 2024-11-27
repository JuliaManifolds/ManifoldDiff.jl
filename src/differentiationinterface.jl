function ManifoldDiff._derivative(f, p, backend::AbstractADType)
    return DI.derivative(f, backend, p)
end

function ManifoldDiff._derivative!(f, X, t, backend::AbstractADType)
    return DI.derivative!(f, X, backend, t)
end

function ManifoldDiff._gradient(f, p, backend::AbstractADType)
    return DI.gradient(f, backend, p)
end

function ManifoldDiff._gradient!(f, X, t, backend::AbstractADType)
    return DI.gradient!(f, X, backend, t)
end

function ManifoldDiff._jacobian(f, p, backend::AbstractADType)
    return DI.jacobian(f, backend, p)
end

function ManifoldDiff._jacobian!(f, X, p, backend::AbstractADType)
    return DI.jacobian!(f, X, backend, p)
end

function ManifoldDiff._hessian(f, p, backend::AbstractADType)
    return DI.hessian(f, backend, p)
end

# function ManifoldDiff._hessian!(f, X, p, backend::AbstractADType)
#     return DI.hessian!(f, X, backend, p)
# end
