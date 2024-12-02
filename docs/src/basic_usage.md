# Basic usage

You can calculate Riemannian gradient of a function defined in its embedding in multiple ways.
For example, `rb_onb_fd51` corresponds to a finite differencing scheme and `rb_onb_fwdd` calculates gradient using `ForwardDiff.jl`. [`DifferentiationInterface.jl`](https://github.com/JuliaDiff/DifferentiationInterface.jl) is used to select the backend.

```@example
using ManifoldDiff
using DifferentiationInterface
using Manifolds, FiniteDifferences, ForwardDiff

rb_onb_fd51 = ManifoldDiff.TangentDiffBackend(AutoFiniteDifferences(central_fdm(5, 1)))
rb_onb_fwdd = ManifoldDiff.TangentDiffBackend(AutoForwardDiff())

s2 = Sphere(2)

A = [1.0 2.0 5.0; 2.0 -1.0 4.0; 5.0 4.0 0.0]

f(p) = p' * A * p
q = [0.0, 1.0, 0.0]

println(ManifoldDiff.gradient(s2, f, q, rb_onb_fd51))
println(ManifoldDiff.gradient(s2, f, q, rb_onb_fwdd))
```

`ManifoldDiff.jl` handles conversion of a Euclidean gradient in the embedding to the Riemannian one.
