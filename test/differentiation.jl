
using Manifolds
using ManifoldDiff
using Test

using ManifoldDiff:
    default_differential_backend,
    _derivative,
    _derivative!,
    differential,
    differential!,
    gradient,
    gradient!,
    _gradient,
    _gradient!,
    _jacobian,
    _jacobian!,
    set_default_differential_backend!,
    AbstractRiemannianDiffBackend

using ManifoldDiff: ExplicitEmbeddedBackend

import ManifoldDiff: gradient

using ADTypes

struct TestRiemannianBackend <: AbstractRiemannianDiffBackend end
function ManifoldDiff.gradient(::AbstractManifold, f, p, ::TestRiemannianBackend)
    return collect(1.0:length(p))
end

using FiniteDifferences: central_fdm
using LinearAlgebra: Diagonal, dot

@testset "Differentiation backend" begin
    fd51 = AutoFiniteDifferences(central_fdm(5, 1))
    @testset "default_differential_backend" begin
        #ForwardDiff is loaded first utils.
        @test default_differential_backend() isa AutoForwardDiff

        @test length(fd51.fdm.grid) == 5
        # check method order
        @test typeof(fd51.fdm).parameters[2] == 1
        fd71 = AutoFiniteDifferences(central_fdm(7, 1))
        @test set_default_differential_backend!(fd71) == fd71
        @test default_differential_backend() == fd71
    end

    using ForwardDiff: ForwardDiff

    fwd_diff = AutoForwardDiff()
    @testset "ForwardDiff" begin
        @test default_differential_backend() isa AutoFiniteDifferences

        @test set_default_differential_backend!(fwd_diff) == fwd_diff
        @test default_differential_backend() == fwd_diff
        @test set_default_differential_backend!(fd51) isa AutoFiniteDifferences
        @test default_differential_backend() isa AutoFiniteDifferences

        set_default_differential_backend!(fwd_diff)
        @test default_differential_backend() == fwd_diff
        set_default_differential_backend!(fd51)
    end

    using FiniteDiff: FiniteDiff

    finite_diff = AutoFiniteDiff()
    @testset "FiniteDiff" begin
        @test default_differential_backend() isa AutoFiniteDifferences

        @test set_default_differential_backend!(finite_diff) == finite_diff
        @test default_differential_backend() == finite_diff
        @test set_default_differential_backend!(fd51) isa AutoFiniteDifferences
        @test default_differential_backend() isa AutoFiniteDifferences

        set_default_differential_backend!(finite_diff)
        @test default_differential_backend() == finite_diff
        set_default_differential_backend!(fd51)
    end

    using ReverseDiff: ReverseDiff

    reverse_diff = AutoReverseDiff()
    @testset "ReverseDiff" begin
        @test default_differential_backend() isa AutoFiniteDifferences

        @test set_default_differential_backend!(reverse_diff) == reverse_diff
        @test default_differential_backend() == reverse_diff
        @test set_default_differential_backend!(fd51) isa AutoFiniteDifferences
        @test default_differential_backend() isa AutoFiniteDifferences

        set_default_differential_backend!(reverse_diff)
        @test default_differential_backend() == reverse_diff
        set_default_differential_backend!(fd51)
    end

    using Zygote: Zygote
    zygote_diff = AutoZygote()

    @testset "gradient" begin
        set_default_differential_backend!(fd51)
        r2 = Euclidean(2)

        c1(t) = [sin(t), cos(t)]
        f1(x) = x[1] + x[2]^2
        function f1!(y, x)
            y .= x[1] + x[2]^2
            return y
        end
        f2(x) = 3 * x[1] * x[2] + x[2]^3
        @test _jacobian(c1, 0.0) ≈ [1.0; 0.0]

        @testset "Inference" begin
            X = [-1.0, -1.0]
            @test (@inferred _derivative(c1, 0.0, AutoForwardDiff())) ≈ [1.0, 0.0]
            @test (@inferred _derivative!(c1, X, 0.0, AutoForwardDiff())) === X
            @test X ≈ [1.0, 0.0]

            @test (@inferred _derivative(c1, 0.0, finite_diff)) ≈ [1.0, 0.0]
            @test (@inferred _gradient(f1, [1.0, -1.0], finite_diff)) ≈ [1.0, -2.0]
        end

        @testset "$(nameof(typeof(backend)))" for backend in [fd51, fwd_diff, finite_diff]
            set_default_differential_backend!(backend)
            @test _derivative(c1, 0.0) ≈ [1.0, 0.0]
            X = [-1.0, -1.0]
            @test _derivative!(c1, X, 0.0) === X
            @test isapprox(X, [1.0, 0.0])
        end
        @testset "$(nameof(typeof(backend)))" for backend in [
            fd51,
            fwd_diff,
            finite_diff,
            reverse_diff,
            zygote_diff,
        ]
            set_default_differential_backend!(backend)
            X = [-1.0, -1.0]
            @test _gradient(f1, [1.0, -1.0]) ≈ [1.0, -2.0]
            @test _gradient!(f1, X, [1.0, -1.0]) === X
            @test X ≈ [1.0, -2.0]
        end
        @testset "$(nameof(typeof(backend)))" for backend in [finite_diff]
            set_default_differential_backend!(backend)
            X = [-0.0 -0.0]
            @test_broken _jacobian(f1, [1.0, -1.0]) ≈ [1.0 -2.0]  # scalar output??
            # The following seems not to work for :central, but it does for forward
            fdf = AutoFiniteDiff()
            @test_broken _jacobian!(f1!, X, [1.0, -1.0], fdf) === X
            @test_broken X ≈ [1.0 -2.0]
        end
        set_default_differential_backend!(ManifoldDiff.NoneDiffBackend())

        @test _jacobian(c1, 0.0, fd51) ≈ [1.0; 0.0]
        jac = [NaN; NaN]
        _jacobian!(c1, jac, 0.0, fd51)
        @test jac ≈ [1.0; 0.0]

        @testset "$(nameof(typeof(backend)))" for backend in [fd51, AutoForwardDiff()]
            @test _derivative(c1, 0.0, backend) ≈ [1.0, 0.0]
            @test _gradient(f1, [1.0, -1.0], backend) ≈ [1.0, -2.0]
        end

        set_default_differential_backend!(fd51)
    end
end

rb_onb_default = TangentDiffBackend(
    default_differential_backend(),
    Manifolds.ExponentialRetraction(),
    Manifolds.LogarithmicInverseRetraction(),
    DefaultOrthonormalBasis(),
    DefaultOrthonormalBasis(),
)

rb_onb_fd51 = TangentDiffBackend(AutoFiniteDifferences(central_fdm(5, 1)))

rb_onb_fwd_diff = TangentDiffBackend(AutoForwardDiff())

rb_onb_finite_diff = TangentDiffBackend(AutoFiniteDiff())

rb_onb_default2 = TangentDiffBackend(
    default_differential_backend();
    basis_arg = CachedBasis(
        DefaultOrthonormalBasis(),
        [[0.0, -1.0, 0.0], [sqrt(2) / 2, 0.0, -sqrt(2) / 2]],
    ),
)

rb_proj = ManifoldDiff.RiemannianProjectionBackend(default_differential_backend())

@testset "Riemannian differentials" begin
    s2 = Sphere(2)
    p = [0.0, 0.0, 1.0]
    q = [1.0, 0.0, 0.0]
    c1(t) = geodesic(s2, q, p, t)

    Xval = [-sqrt(2) / 2, 0.0, sqrt(2) / 2]
    @test isapprox(s2, c1(π / 4), differential(s2, c1, π / 4, rb_onb_default), Xval)
    X = similar(p)
    differential!(s2, c1, X, π / 4, rb_onb_default)
    @test isapprox(s2, c1(π / 4), X, Xval)

    @testset for backend in [rb_onb_fd51, rb_onb_fwd_diff, rb_onb_finite_diff]
        @test isapprox(s2, c1(π / 4), differential(s2, c1, π / 4, backend), Xval)
        X = similar(p)
        differential!(s2, c1, X, π / 4, backend)
        @test isapprox(s2, c1(π / 4), X, Xval)
    end
end

@testset "Riemannian gradients" begin
    s2 = Sphere(2)
    f1(p) = p[1]

    q = [sqrt(2) / 2, 0, sqrt(2) / 2]
    X = similar(q)
    for backend in [rb_onb_default, rb_onb_default2, rb_proj]
        @test isapprox(s2, q, gradient(s2, f1, q, backend), [0.5, 0.0, -0.5])
        @test gradient!(s2, f1, X, q, backend) === X
        @test isapprox(s2, q, X, [0.5, 0.0, -0.5])
    end
    X = similar(q)
    for backend in [rb_onb_default, rb_onb_default2, rb_proj]
        gradient!(s2, f1, X, q, backend)
        @test isapprox(s2, q, X, [0.5, 0.0, -0.5])
    end

    # Test the gradient fallback
    @test gradient(s2, f1, q, TestRiemannianBackend()) == [1.0, 2.0, 3.0]
    X = similar(q)
    @test gradient!(s2, f1, X, q, TestRiemannianBackend()) === X
    @test X == [1.0, 2.0, 3.0]
end

@testset "Riemannian Jacobians" begin
    s2 = Sphere(2)
    f1(p) = p

    q = [sqrt(2) / 2, 0, sqrt(2) / 2]
    X = similar(q)
    @test isapprox(
        s2,
        q,
        ManifoldDiff.jacobian(s2, s2, f1, q, rb_onb_default),
        [1.0 0.0; 0.0 1.0],
    )

    q2 = [1.0, 0.0, 0.0]
    f2(X) = [0.0 0.0 0.0; 0.0 2.0 -1.0; 0.0 -3.0 1.0] * X
    Tq2s2 = TangentSpace(s2, q2)
    @test isapprox(
        ManifoldDiff.jacobian(Tq2s2, Tq2s2, f2, zero_vector(s2, q2), rb_onb_default),
        [2.0 -1.0; -3.0 1.0],
    )

    q3 = [0.0, 1.0, 0.0]
    f3(X) = [0.0 2.0 1.0; 0.0 0.0 0.0; 0.0 5.0 1.0] * X
    Tq3s2 = TangentSpace(s2, q3)
    @test isapprox(
        ManifoldDiff.jacobian(Tq2s2, Tq3s2, f3, zero_vector(s2, q2), rb_onb_default),
        [-2.0 -1.0; 5.0 1.0],
    )
end

@testset "Riemannian Hessians" begin
    s2 = Sphere(2)
    q = [sqrt(2) / 2, 0, sqrt(2) / 2]
    q2 = [0.0, 1.0, 0.0]

    f1(p) = distance(s2, q2, p)^2

    @test isapprox(ManifoldDiff.hessian(s2, f1, q, rb_onb_default), [2.0 0.0; 0.0 0.0])
    @test_broken isapprox(ManifoldDiff.hessian(s2, f1, q, rb_onb_fwd_diff), [2.0 0.0; 0.0 0.0])
end

@testset "EmbeddedBackend" begin
    A = [1 0 0; 0 2 0; 0 0 3.0]
    p = 1 / sqrt(2.0) .* [1.0, 1.0, 0.0]

    cost = p -> p' * A * p # Euclidean cost of Rayleigh
    grad = (M, p) -> 2A * p
    grad! = (M, X, p) -> X .= 2A * p
    M = Euclidean(3)
    E = ExplicitEmbeddedBackend(M, gradient = grad, (gradient!) = grad!)
    S = Sphere(2)
    r_grad = (S, p) -> project(S, p, grad(M, p))
    Xt = r_grad(S, p)
    @test is_vector(S, p, Xt, true; atol = 1e-14)

    R = RiemannianProjectionBackend(E)
    X = gradient(S, cost, p, R)
    @test isapprox(S, p, X, Xt)
    gradient!(S, cost, X, p, R)
    @test isapprox(S, p, X, Xt)
    # Errors with empty gradient
    Ee = ExplicitEmbeddedBackend(M)
    Re = RiemannianProjectionBackend(Ee)
    @test_throws MissingException gradient(S, cost, p, Re)
    @test_throws MissingException gradient!(S, cost, X, p, Re)
end
