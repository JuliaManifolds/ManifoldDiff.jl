module ManifoldDiffFiniteDifferencesExt

using ManifoldDiff
using ManifoldDiff: FiniteDifferencesBackend
using FiniteDifferences

function ManifoldDiff.FiniteDifferencesBackend()
    return FiniteDifferencesBackend(central_fdm(5, 1))
end

function ManifoldDiff._derivative(f, t, backend::FiniteDifferencesBackend)
    return backend.method(f, t)
end

function ManifoldDiff._gradient(f, p, backend::FiniteDifferencesBackend)
    return FiniteDifferences.grad(backend.method, f, p)[1]
end

function ManifoldDiff._jacobian(f, p, backend::FiniteDifferencesBackend)
    return FiniteDifferences.jacobian(backend.method, f, p)[1]
end

function ManifoldDiff._hessian(f, p, backend::FiniteDifferencesBackend)
    return FiniteDifferences.jacobian(
        backend.method,
        q -> FiniteDifferences.grad(backend.method, f, q)[1],
        p,
    )[1]
end

end
