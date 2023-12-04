using Manifolds, ManifoldDiff, Test

@testset "proximal maps" begin
    #
    # Distance
    p = [1.0, 0.0, 0.0]
    q = [0.0, 1.0, 0.0]
    M = Sphere(2)
    N = PowerManifold(M, NestedPowerRepresentation(), 2)
    @test_throws ErrorException ManifoldDiff.prox_distance(M, 1.0, p, q, 3)
    @test_throws ErrorException ManifoldDiff.prox_distance!(M, p, 1.0, p, q, 3)
    @test distance(
        M,
        ManifoldDiff.prox_distance(M, distance(M, p, q) / 2, p, q, 1),
        shortest_geodesic(M, p, q, 0.5),
    ) < eps()
    t = similar(p)
    ManifoldDiff.prox_distance!(M, t, distance(M, p, q) / 2, p, q, 1)
    @test t == ManifoldDiff.prox_distance(M, distance(M, p, q) / 2, p, q, 1)
end
