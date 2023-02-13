using Manifolds, ManifoldsBase, ManifoldDiff, Test

@testset "Derivatives" begin
    M = Sphere(2)
    p = [0.0, 0.0, 1.0]
    q = [0.0, 1.0, 0.0]
    r = [1.0, 0.0, 0.0]
    @testset "Geodesic" begin
        X = log(M, p, q)
        s = geodesic(M, p, X, 0.5)
        Y = zero_vector(M, p)
        ManifoldDiff.geodesic_derivative!(M, Y, p, X, 0.5)
        Z = ManifoldDiff.geodesic_derivative(M, p, X, 0.5)
        sol = parallel_transport_to(M, p, X, s)
        @test isapprox(M, s, sol, Y)
        @test isapprox(M, s, sol, Z)

        ManifoldDiff.shortest_geodesic_derivative!(M, Y, p, q, 0.5)
        Z = ManifoldDiff.shortest_geodesic_derivative(M, p, q, 0.5)
        sol = parallel_transport_to(M, p, X, s)
        @test isapprox(M, s, sol, Y)
        @test isapprox(M, s, sol, Z)
    end
end
