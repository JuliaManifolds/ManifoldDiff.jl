using ManifoldDiff
using Manifolds
using ManifoldsBase
using LinearAlgebra
using Test

@testset "ManifoldDiff GPU AD gradients" begin
    cuda_loaded = false
    try
        using CUDA
        cuda_loaded = CUDA.functional()
    catch
        cuda_loaded = false
    end

    if cuda_loaded
        @eval using CUDA

        # Note: ForwardDiff does NOT support CuArrays (scalar indexing in seed!).
        # Only Zygote (reverse-mode) is tested here.
        zygote_loaded = false
        try
            using Zygote
            using ADTypes: AutoZygote
            zygote_loaded = true
        catch
            zygote_loaded = false
        end

        if zygote_loaded
            @eval using Zygote
            @eval using ADTypes: AutoZygote

            @testset "Euclidean gradient — Zygote Float64" begin
                M = Euclidean(3)
                backend = ManifoldDiff.RiemannianProjectionBackend(AutoZygote())

                # f(p) = sum(p.^2) / 2, analytical grad = p
                f(p) = sum(p .^ 2) / 2

                p_cpu = Float64[1.0, 2.0, 3.0]
                p = CuArray(p_cpu)

                grad = ManifoldDiff.gradient(M, f, p, backend)
                @test grad isa CuArray{Float64}
                @test isapprox(Array(grad), p_cpu; atol=1e-10)
            end

            @testset "Euclidean gradient — Zygote Float32" begin
                M = Euclidean(3)
                backend = ManifoldDiff.RiemannianProjectionBackend(AutoZygote())

                f(p) = sum(p .^ 2) / 2f0

                p_cpu = Float32[1.0, 2.0, 3.0]
                p = CuArray(p_cpu)

                grad = ManifoldDiff.gradient(M, f, p, backend)
                @test grad isa CuArray{Float32}
                @test isapprox(Array(grad), p_cpu; atol=1e-4)
            end

            @testset "Sphere gradient — Zygote" begin
                M = Sphere(2)
                backend = ManifoldDiff.RiemannianProjectionBackend(AutoZygote())

                # Use dot product instead of p[1] to avoid scalar indexing.
                # f(p) = dot(a, p) where a = [1,0,0]
                # Euclidean grad = a, Riemannian grad = a - dot(a,p)*p
                a = CuArray([1.0, 0.0, 0.0])
                f(p) = dot(a, p)

                p_cpu = [sqrt(2.0) / 2, 0.0, sqrt(2.0) / 2]
                p = CuArray(p_cpu)
                a_cpu = [1.0, 0.0, 0.0]
                expected_grad = a_cpu - dot(a_cpu, p_cpu) * p_cpu

                grad = ManifoldDiff.gradient(M, f, p, backend)
                @test grad isa CuArray{Float64}
                # Should be in tangent space: dot(grad, p) ≈ 0
                @test abs(dot(Array(grad), p_cpu)) < 1e-10
                @test isapprox(Array(grad), expected_grad; atol=1e-10)
            end

            @testset "CPU vs GPU gradient equivalence — Zygote" begin
                M = Euclidean(5)
                backend = ManifoldDiff.RiemannianProjectionBackend(AutoZygote())

                f(p) = sum(p .^ 2) / 2
                p_cpu = randn(5)
                p_gpu = CuArray(p_cpu)

                grad_cpu = ManifoldDiff.gradient(M, f, p_cpu, backend)
                grad_gpu = ManifoldDiff.gradient(M, f, p_gpu, backend)
                @test grad_gpu isa CuArray{Float64}
                @test isapprox(Array(grad_gpu), grad_cpu; atol=1e-12)
            end

            @testset "Quadratic objective — Zygote" begin
                M = Euclidean(4)
                backend = ManifoldDiff.RiemannianProjectionBackend(AutoZygote())

                target = CuArray([1.0, 2.0, 3.0, 4.0])
                f(p) = sum((p .- target) .^ 2) / 2

                p_cpu = zeros(4)
                p = CuArray(p_cpu)

                grad = ManifoldDiff.gradient(M, f, p, backend)
                @test grad isa CuArray{Float64}
                # grad = p - target = [0,0,0,0] - [1,2,3,4] = [-1,-2,-3,-4]
                @test isapprox(Array(grad), p_cpu .- Array(target); atol=1e-10)
            end
        else
            @info "Zygote not available, skipping GPU AD tests"
        end
    else
        @info "CUDA not functional, skipping ManifoldDiff GPU AD tests"
    end
end
