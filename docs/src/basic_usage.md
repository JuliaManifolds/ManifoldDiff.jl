# Basic usage

You can calculate Riemannian gradient of a function defined in its embedding in multiple ways.
[`DifferentiationInterface.jl`](https://github.com/JuliaDiff/DifferentiationInterface.jl) can be used to select the backend.

```@example
using ManifoldDiff
using DifferentiationInterface
using Manifolds, FiniteDifferences, ForwardDiff, Zygote

rb_onb_fd51 = TangentDiffBackend(AutoFiniteDifferences(central_fdm(5, 1)))
rb_onb_fwdd = TangentDiffBackend(AutoForwardDiff())
rb_proj_zyg = RiemannianProjectionBackend(AutoZygote())

s2 = Sphere(2)

A = [1.0 2.0 5.0; 2.0 -1.0 4.0; 5.0 4.0 0.0]

f(p) = p' * A * p
q = [0.0, 1.0, 0.0]

println(ManifoldDiff.gradient(s2, f, q, rb_onb_fd51))
println(ManifoldDiff.gradient(s2, f, q, rb_onb_fwdd))
println(ManifoldDiff.gradient(s2, f, q, rb_proj_zyg))
```

In this example `rb_onb_fd51` corresponds to a finite differencing scheme, `rb_onb_fwdd` calculates gradient using [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl) and `rb_proj_zyg` uses [`Zygote.jl`](https://github.com/FluxML/Zygote.jl) for reverse mode automatic differentiation.

[`TangentDiffBackend`](@ref) reduces dimensionality of the problem to the intrinsic dimension of the manifold, while [`RiemannianProjectionBackend`](@ref) relies on converting Euclidean gradient in the embedding to the Riemannian one.
