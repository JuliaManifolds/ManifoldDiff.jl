using Manifolds, ManifoldsBase, ManifoldDiff, Test
using FiniteDifferences

function finitedifferences_jacobian(
    M_arg::AbstractManifold,
    M_val::AbstractManifold,
    f,
    p;
    B_arg::AbstractBasis = DefaultOrthonormalBasis(),
    B_val::AbstractBasis = DefaultOrthonormalBasis(),
    m = central_fdm(5, 1),
)
    fv = f(p)
    function fj(_c)
        Y = get_vector(M_arg, p, _c, B_arg)
        return get_coordinates(M_val, fv, log(M_val, fv, f(exp(M_arg, p, Y))), B_val)
    end
    return FiniteDifferences.jacobian(m, fj, zeros(manifold_dimension(M_arg)))[1]
end

@testset "Jacobians" begin
    M = Sphere(2)
    p = [0.0, 0.0, 1.0]
    X = [0.0, 1.0, 0.0]
    q = [0.0, 1.0, 0.0]

    Jfnd = finitedifferences_jacobian(TangentSpace(M, p), M, Y -> exp(M, p, Y), X)

    J = ManifoldDiff.jacobian_exp_argument(M, p, X)
    @test J ≈ Jfnd atol = 1e-8

    ManifoldDiff.differential_exp_argument(M, p, X, [1.0, 0.0, 0.0])

    Jc = similar(J)
    ManifoldDiff.jacobian_exp_argument!(M, Jc, p, X)
    @test Jc ≈ Jfnd atol = 1e-8

    Jfnd = finitedifferences_jacobian(
        M,
        M,
        q -> exp(M, q, parallel_transport_to(M, p, X, q)),
        p,
    )
    J = ManifoldDiff.jacobian_exp_basepoint(M, p, X)
    @test J ≈ Jfnd atol = 1e-8

    Jc = similar(J)
    ManifoldDiff.jacobian_exp_basepoint!(M, Jc, p, X)
    @test Jc ≈ Jfnd atol = 1e-8

    Jfnd = finitedifferences_jacobian(M, TangentSpace(M, p), q -> log(M, p, q), q)

    J = ManifoldDiff.jacobian_log_argument(M, p, q)
    @test J ≈ Jfnd atol = 1e-8

    Jc = similar(J)
    ManifoldDiff.jacobian_log_argument!(M, Jc, p, q)
    @test Jc ≈ Jfnd atol = 1e-8

    Jref = finitedifferences_jacobian(
        M,
        TangentSpace(M, p),
        r -> parallel_transport_to(M, r, log(M, r, q), p),
        p,
    )
    J = ManifoldDiff.jacobian_log_basepoint(M, p, q)
    @test J ≈ Jref atol = 1e-8

    Jc = similar(J)
    ManifoldDiff.jacobian_log_basepoint!(M, Jc, p, q)
    @test Jc ≈ Jref atol = 1e-8
end
