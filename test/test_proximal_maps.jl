using Manifolds, ManifoldDiff, Test

@testset "proximal maps" begin
    #
    # Distance
    p = [1.0, 0.0, 0.0]
    q = [0.0, 1.0, 0.0]
    M = Sphere(2)
    # Error for r != 1 or 2
    @test_throws ErrorException ManifoldDiff.prox_distance(M, 1.0, p, q, 3)
    @test_throws ErrorException ManifoldDiff.prox_distance!(M, p, 1.0, p, q, 3)
    # r = 1
    @test distance(
        M,
        ManifoldDiff.prox_distance(M, distance(M, p, q) / 2, p, q, 1),
        shortest_geodesic(M, p, q, 0.5),
    ) < eps()
    t = similar(p)
    ManifoldDiff.prox_distance!(M, t, distance(M, p, q) / 2, p, q, 1)
    @test t == ManifoldDiff.prox_distance(M, distance(M, p, q) / 2, p, q, 1)
    # r = 2
    ManifoldDiff.prox_distance!(M, t, 1.0, p, q, 2)
    @test t == ManifoldDiff.prox_distance(M, 1.0, p, q, 2)
    @test t == ManifoldDiff.prox_distance(M, 1.0, p, q) # 2 is also the default
    @test distance(M, t, shortest_geodesic(M, p, q, 0.5)) ≈ 0 atol = 1e-15
end
