# Differentiation backends

```@docs
set_default_differential_backend!
default_differential_backend
```

## Euclidian backends

Euclidian backend objects can be taken from [ADTypes.jl](https://github.com/SciML/ADTypes.jl).
See the documentation of [DifferentiationInterface.jl](https://github.com/JuliaDiff/DifferentiationInterface.jl) for the list of supported packages.

## EmbeddedDiff

```@autodocs
Modules = [ManifoldDiff]
Pages = ["embedded_diff.jl"]
Order = [:type, :function, :constant]
```
