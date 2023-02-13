
using Manifolds, ManifoldDiff
using Test
using LinearAlgebra

using ManifoldDiff:
    differential_exp_argument,
    differential_exp_argument_lie_approx,
    differential_inverse_retract_argument_fd_approx

@testset "Rotations(3)" begin
    M = Rotations(3)
    p = rand(M)
    X = rand(M; vector_at = p)
    Y = rand(M; vector_at = p)

    @test differential_exp_argument(M, p, X, Y) ≈
          invoke(differential_exp_argument, Tuple{AbstractManifold,Any,Any,Any}, M, p, X, Y)


    G = SpecialOrthogonal(3)
    p = Matrix(I, 3, 3)

    ω = [[1.0, 2.0, 3.0], [3.0, 2.0, 1.0], [1.0, 3.0, 2.0]]
    pts = [exp(M, p, hat(M, p, ωi)) for ωi in ω]
    Xpts = [hat(M, p, [-1.0, 2.0, 0.5]), hat(M, p, [1.0, 0.0, 0.5])]

    @testset "differentials" begin
        q2 = exp(G, pts[1], Xpts[2])
        @test isapprox(
            G,
            q2,
            differential_exp_argument_lie_approx(G, pts[1], Xpts[1], Xpts[2]; n = 0),
            Xpts[2],
        )
        diff_ref = [
            0.0 -0.7482721017619345 -0.508151233069837
            0.7482721017619345 0.0 -0.10783358474129323
            0.508151233069837 0.10783358474129323 0.0
        ]
        @test isapprox(
            G,
            q2,
            differential_exp_argument_lie_approx(G, pts[1], Xpts[1], Xpts[2]),
            diff_ref;
            atol = 1e-12,
        )
    end
end

@testset "FiniteDifferenceLogDiffArgumentMethod" begin
    M = Sphere(2)
    p = [1.0, 0.0, 0.0]
    q = [0.0, sqrt(2) / 2, sqrt(2) / 2]
    X = [1.0, -2.0, 2.0]

    # computed using Manopt.differential_log_argument(M, p, q, X)
    diff_ref = [-5.131524956784507e-33, -3.84869943477634, 2.434485872403245]
    @test isapprox(
        M,
        q,
        differential_inverse_retract_argument_fd_approx(M, p, q, X; h = 1e-4),
        diff_ref;
        atol = 1e-7,
    )
end
