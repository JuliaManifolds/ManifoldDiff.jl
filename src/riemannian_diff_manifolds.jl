
@doc raw"""
    riemannian_gradient(M, p, Y; embedding_metric=EuclideanMetric())
    riemannian_gradient!(M, X, p, Y; embedding_metric=EuclideanMetric())

For a given gradient ``Y = \operatorname{grad} \tilde f(p)`` in the embedding of a manifold,
this function computes the Riemannian gradient ``\operatorname{grad} f(p)`` of the function
``\tilde f`` restricted to the manifold ``M``.
This can also be done in place of `X`.

By default it uses the following computation:
Let the projection ``Z = \operatorname{proj}_{T_p\mathcal M}(Y)`` of ``Y`` onto the
tangent space at ``p`` be given, that is with respect to the inner product in the embedding.
Then

```math
⟨Z-Y, W⟩ = 0 \text{ for all } W \in T_p\mathcal M,
```

or rearranged ``⟨Y,W⟩ = ⟨Z,W⟩``. We then use the definition of the Riemannian gradient

```math
⟨\operatorname{grad} f(p), W⟩_p = Df(p)[X] = ⟨\operatorname{grad}f(p), W⟩ = ⟨\operatorname{proj}_{T_p\mathcal M}(\operatorname{grad}f(p)),W⟩
\quad\text{for all } W \in T_p\mathcal M.
```
Comparing the first and the last term, the remaining computation is the function [`change_representer`](@ref change_representer(M::AbstractManifold, G2::AbstractMetric, p, X)).

This method can also be implemented directly, if a more efficient/stable version is known.

The function is inspired by `egrad2rgrad` in the [Matlab package Manopt](https://manopt.org).
"""
function riemannian_gradient(
    M::AbstractManifold,
    p,
    Y;
    embedding_metric::AbstractMetric = EuclideanMetric(),
)
    X = zero_vector(M, p)
    riemannian_gradient!(M, X, p, Y; embedding_metric = embedding_metric)
    return X
end

function riemannian_gradient!(
    M::AbstractManifold,
    X,
    p,
    Y;
    embedding_metric = EuclideanMetric(),
)
    project!(M, X, p, Y)
    change_representer!(M, X, embedding_metric, p, X)
    return X
end
