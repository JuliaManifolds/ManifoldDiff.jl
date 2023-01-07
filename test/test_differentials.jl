using Manifolds, ManifoldDiff
using Test

using ManifoldDiff:
    differential_exp_argument,
    differential_exp_argument!,
    differential_exp_basepoint,
    differential_exp_basepoint!,
    differential_log_argument,
    differential_log_argument!,
    differential_log_basepoint,
    differential_log_basepoint!,
    differential_shortest_geodesic_endpoint,
    differential_shortest_geodesic_endpoint!,
    differential_shortest_geodesic_startpoint,
    differential_shortest_geodesic_startpoint!


@testset "Differentials" begin
    p = [1.0, 0.0, 0.0]
    q = [0.0, 1.0, 0.0]
    M = Sphere(2)
    X = log(M, p, q)
    X_qp = log(M, q, p)
    Y = similar(X)
    @testset "Differentials on Sn(2)" begin
        @test differential_log_basepoint(M, p, p, X) == -X
        differential_log_basepoint!(M, Y, p, p, X)
        @test Y == -X
        @test differential_log_basepoint(M, p, q, X) == -X
        differential_log_basepoint!(M, Y, p, q, X) == -X
        @test Y == -X
        @test differential_log_argument(M, p, p, X) == X
        differential_log_argument!(M, Y, p, p, X) == X
        @test Y == X
        @test differential_log_argument(M, p, q, X_qp) == -X
        differential_log_argument!(M, Y, p, q, X_qp)
        @test Y == -X
        @test differential_exp_basepoint(M, p, zero_vector(M, p), X) == X
        differential_exp_basepoint!(M, Y, p, zero_vector(M, p), X)
        @test Y == X
        @test norm(M, q, differential_exp_basepoint(M, p, X, X) - [-π / 2, 0.0, 0.0]) ≈ 0 atol =
            6 * 10^(-16)
        differential_exp_basepoint!(M, Y, p, X, X)
        @test norm(M, q, Y - [-π / 2, 0.0, 0.0]) ≈ 0 atol = 6 * 10^(-16)
        @test differential_exp_argument(M, p, zero_vector(M, p), X) == X
        differential_exp_argument!(M, Y, p, zero_vector(M, p), X) == X
        @test Y == X
        @test norm(M, q, differential_exp_argument(M, p, X, zero_vector(M, p))) ≈ 0
        differential_exp_argument!(M, Y, p, X, zero_vector(M, p))
        @test norm(M, q, Y) ≈ 0
        for t in [0, 0.15, 0.33, 0.66, 0.9]
            @test differential_shortest_geodesic_startpoint(M, p, p, t, X) == (1 - t) * X
            differential_shortest_geodesic_startpoint!(M, Y, p, p, t, X)
            @test Y == (1 - t) * X
            @test norm(
                M,
                p,
                differential_shortest_geodesic_endpoint(M, p, p, t, X) - t * X,
            ) ≈ 0 atol = 10.0^(-16)
            differential_shortest_geodesic_endpoint!(M, Y, p, p, t, X)
            @test norm(M, p, Y - t * X) ≈ 0 atol = 10.0^(-16)
        end
    end
    @testset "Differentials on Power of Sn(2)" begin
        for rep in [NestedPowerRepresentation(), NestedReplacingPowerRepresentation()]
            N = PowerManifold(M, rep, 3)
            pN = [p, q, p]
            qN = [p, p, q]
            V = [X, zero_vector(M, p), -X_qp]
            W = allocate(V)
            @test differential_log_argument(N, pN, qN, V) == [V[1], V[2], V[1]]
            differential_log_argument!(N, W, pN, qN, V)
            @test W == [V[1], V[2], V[1]]
        end
    end
    @testset "Differentials on SPD(2)" begin
        #
        # Single differentials on SPD
        M2 = SymmetricPositiveDefinite(2)
        p2 = [1.0 0.0; 0.0 1.0]
        X2 = [0.5 1.0; 1.0 0.5]
        q2 = exp(M2, p2, X2)
        # Test differentials (1) Dx of Log_xy
        @test norm(M2, p2, differential_log_basepoint(M2, p2, p2, X2) + X2) ≈ 0 atol =
            4 * 10^(-16)
        @test norm(M2, q2, differential_log_argument(M2, p2, q2, zero_vector(M2, p2))) ≈ 0 atol =
            4 * 10^(-16)
        @test norm(
            M2,
            p2,
            differential_exp_basepoint(M2, p2, zero_vector(M2, p2), X2) - X2,
        ) ≈ 0 atol = 4 * 10^(-16)
        @test norm(
            M2,
            p2,
            differential_exp_argument(M2, p2, zero_vector(M2, p2), X2) - X2,
        ) ≈ 0 atol = 4 * 10^(-16)

        @test norm(M2, q2, differential_exp_basepoint(M2, p2, X2, zero_vector(M2, p2))) ≈ 0 atol =
            4 * 10.0^(-16)
        @test norm(M2, q2, differential_exp_argument(M2, p2, X2, zero_vector(M2, p2))) ≈ 0 atol =
            4 * 10.0^(-16)
        # test coeff of log_basepoint, since it is not always expicitly used.
        @test ManifoldDiff.βdifferential_log_basepoint(-1.0, 1.0, 2.0) ≈
              -2 * cosh(2.0) / sinh(2.0)
    end
    @testset "Differentials on Euclidean(2)" begin
        M3 = Euclidean(2)
        x3 = [1.0, 2.0]
        ξ3 = [1.0, 0.0]
        @test norm(M3, x3, differential_exp_basepoint(M3, x3, ξ3, ξ3) - ξ3) ≈ 0 atol =
            4 * 10.0^(-16)
    end
    @testset "Differentials on the Circle" begin
        M = Circle()
        p = 0
        q = π / 4
        X = π / 8
        @test differential_log_argument(M, p, q, X) == X
    end
    @testset "Differentials on Euclidean()" begin
        M = Euclidean()
        p = 0
        q = π / 4
        X = π / 8
        @test differential_log_argument(M, p, q, X) == X
    end

    @testset "Differentials on ProductManifold" begin
        M = ProductManifold(Sphere(2), SymmetricPositiveDefinite(2))
        p = ArrayPartition([0.0, 0.0, 1.0], [1.0 0.0; 0.0 1.0])
        Xp = ArrayPartition([-0.5, 1.0, 0.0], [0.5 1.0; 1.0 0.5])
        q = exp(M, p, Xp)
        # Test differentials (1) Dx of Log_xy
        @test norm(M, p, differential_log_basepoint(M, p, p, Xp) + Xp) ≈ 0 atol =
            4 * 10^(-16)
        @test norm(M, q, differential_log_argument(M, p, q, zero_vector(M, p))) ≈ 0 atol =
            4 * 10^(-16)
        @test norm(M, p, differential_exp_basepoint(M, p, zero_vector(M, p), Xp) - Xp) ≈ 0 atol =
            4 * 10^(-16)
        @test norm(M, p, differential_exp_argument(M, p, zero_vector(M, p), Xp) - Xp) ≈ 0 atol =
            4 * 10^(-16)

        @test norm(M, q, differential_exp_basepoint(M, p, Xp, zero_vector(M, p))) ≈ 0 atol =
            4 * 10.0^(-16)
        @test norm(M, q, differential_exp_argument(M, p, Xp, zero_vector(M, p))) ≈ 0 atol =
            4 * 10.0^(-16)
    end
end
