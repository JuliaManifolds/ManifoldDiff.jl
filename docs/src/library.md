# Different library functions

Documentation for `ManifoldDiff.jl`'s methods and types for finite differences and automatic differentiation.

## Derivatives

```@autodocs
Modules = [ManifoldDiff]
Pages = ["derivatives.jl"]
Order = [:type, :function, :constant]
Private = true
```

## Differentials and their adjoints

```@autodocs
Modules = [ManifoldDiff]
Pages = ["adjoint_differentials.jl","differentials.jl"]
Order = [:type, :function, :constant]
Private = true
```

```@autodocs
Modules = [ManifoldDiff]
Pages = ["diagonalizing_projectors.jl"]
Order = [:type, :function, :constant]
Private = true
```

## Gradients

```@autodocs
Modules = [ManifoldDiff]
Pages = ["gradients.jl"]
Order = [:type, :function, :constant]
Private = true
```

## Jacobi fields

```@autodocs
Modules = [ManifoldDiff]
Pages = ["Jacobi_fields.jl"]
Order = [:type, :function, :constant]
```

## Jacobians

```@autodocs
Modules = [ManifoldDiff]
Pages = ["jacobians.jl"]
Order = [:type, :function, :constant]
```

## Riemannian differentials

```@autodocs
Modules = [ManifoldDiff]
Pages = ["riemannian_diff.jl"]
Order = [:type, :function, :constant]
```

## Proximal Maps

Given a convex, lower semi-continuous function ``f\colon \mathcal M \to \mathbb R``, its proximal map is defined
for some ``λ>0`` as [Bacak:2014](@cite)

```math
\operatorname{prox}_{λf}(p) := \operatorname*{arg\,min}_{q\in\mathcal M} \frac{1}{2λ}d^2_{\mathcal M}(p,q) + f(q).
```

Another name for the proximal map is _resolvent_.
Intuitively this means to minimize the function ``f`` while at the same timme “staying close”
to the argument ``p``.

```@autodocs
Modules = [ManifoldDiff]
Pages = ["proximal_maps.jl"]
Order = [:type, :function, :constant]
```
