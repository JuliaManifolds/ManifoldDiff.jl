var documenterSearchIndex = {"docs":
[{"location":"library.html#Different-library-functions","page":"Library","title":"Different library functions","text":"","category":"section"},{"location":"library.html","page":"Library","title":"Library","text":"Documentation for ManifoldDiff.jl's methods and types for finite differences and automatic differentiation.","category":"page"},{"location":"library.html#Derivatives","page":"Library","title":"Derivatives","text":"","category":"section"},{"location":"library.html","page":"Library","title":"Library","text":"Modules = [ManifoldDiff]\nPages = [\"derivatives.jl\"]\nOrder = [:type, :function, :constant]\nPrivate = true","category":"page"},{"location":"library.html#ManifoldDiff.geodesic_derivative-Tuple{Any, Any, Any, Number}","page":"Library","title":"ManifoldDiff.geodesic_derivative","text":"Y = geodesic_derivative(M, p, X, t::Number; γt = geodesic(M, p, X, t))\ngeodesic_derivative!(M, Y, p, X, t::Number; γt = geodesic(M, p, X, t))\n\nEvaluate the derivative of the geodesic γ(t) with γ_pX(0) = p and dot γ_pX(0) = X at t. The formula reads\n\ndot γ(t) = mathcal P_γ(t) gets p X\n\nwhere mathcal P denotes the parallel transport. This computation can also be done in-place of Y.\n\nOptional Parameters\n\nγt – (geodesic(M, p, X, t)) the point on the geodesic at t. This way if the point was computed earlier it can be resued here.\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.shortest_geodesic_derivative-Tuple{Any, Any, Any, Number}","page":"Library","title":"ManifoldDiff.shortest_geodesic_derivative","text":"Y = shortest_geodesic_derivative(M, p, X, t::Number; γt = shortest_geodesic(M, p, q, t))\nshortest_geodesic_derivative!(M, Y, p, X, t::Number; γt = shortest_geodesic(M, p, q, t))\n\nEvaluate the derivative of the shortest geodesic γ(t) with γ_pq(0) = p and dot γ_pq(1) = q at t. The formula reads\n\ndot γ(t) = mathcal P_γ(t) gets p log_pq\n\nwhere mathcal P denotes the parallel transport. This computation can also be done in-place of Y.\n\nOptional Parameters\n\nγt = geodesic(M, p, X, t) the point on the geodesic at t. This way if the point was computed earlier it can be resued here.\n\n\n\n\n\n","category":"method"},{"location":"library.html#Differentials-and-their-adjoints","page":"Library","title":"Differentials and their adjoints","text":"","category":"section"},{"location":"library.html","page":"Library","title":"Library","text":"Modules = [ManifoldDiff]\nPages = [\"adjoint_differentials.jl\"]\nOrder = [:type, :function, :constant]\nPrivate = true","category":"page"},{"location":"library.html#ManifoldDiff.adjoint_differential_exp_argument-Tuple{AbstractManifold, Any, Any, Any}","page":"Library","title":"ManifoldDiff.adjoint_differential_exp_argument","text":"adjoint_differential_exp_argument(M, p, X, Y)\nadjoint_differential_exp_argument!(M, Z, p, X, Y)\n\nCompute the adjoint of D_Xexp_p XY (in place of Z). Note that X   T_p(T_pmathcal M) = T_pmathcal M is still a tangent vector.\n\nSee also\n\ndifferential_exp_argument, adjoint_Jacobi_field\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.adjoint_differential_exp_basepoint-Tuple{AbstractManifold, Any, Any, Any}","page":"Library","title":"ManifoldDiff.adjoint_differential_exp_basepoint","text":"adjoint_differential_exp_basepoint(M, p, X, Y)\nadjoint_differential_exp_basepoint!(M, Z, p, X, Y)\n\nComputes the adjoint of D_p exp_p XY (in place of Z).\n\nSee also\n\ndifferential_exp_basepoint, adjoint_Jacobi_field\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.adjoint_differential_log_argument-Tuple{AbstractManifold, Any, Any, Any}","page":"Library","title":"ManifoldDiff.adjoint_differential_log_argument","text":"adjoint_differential_log_argument(M, p, q, X)\nadjoint_differential_log_argument!(M, Y, p, q, X)\n\nCompute the adjoint of D_q log_p qX (in place of Y).\n\nSee also\n\ndifferential_log_argument, adjoint_Jacobi_field\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.adjoint_differential_log_basepoint-Tuple{AbstractManifold, Any, Any, Any}","page":"Library","title":"ManifoldDiff.adjoint_differential_log_basepoint","text":"adjoint_differential_log_basepoint(M, p, q, X)\nadjoint_differential_log_basepoint!(M, Y, p, q, X)\n\ncomputes the adjoint of D_p log_p qX (in place of Y).\n\nSee also\n\ndifferential_log_basepoint, adjoint_Jacobi_field\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.adjoint_differential_shortest_geodesic_endpoint-Tuple{AbstractManifold, Any, Any, Any, Any}","page":"Library","title":"ManifoldDiff.adjoint_differential_shortest_geodesic_endpoint","text":"adjoint_differential_shortest_geodesic_endpoint(M, p, q, t, X)\nadjoint_differential_shortest_geodesic_endpoint!(M, Y, p, q, t, X)\n\nCompute the adjoint of D_q γ(t p q)X (in place of Y).\n\nSee also\n\ndifferential_shortest_geodesic_endpoint, adjoint_Jacobi_field\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.adjoint_differential_shortest_geodesic_startpoint-Tuple{AbstractManifold, Any, Any, Any, Any}","page":"Library","title":"ManifoldDiff.adjoint_differential_shortest_geodesic_startpoint","text":"adjoint_differential_shortest_geodesic_startpoint(M, p, q, t, X)\nadjoint_differential_shortest_geodesic_startpoint!(M, Y, p, q, t, X)\n\nCompute the adjoint of D_p γ(t p q)X (in place of Y).\n\nSee also\n\ndifferential_shortest_geodesic_startpoint, adjoint_Jacobi_field\n\n\n\n\n\n","category":"method"},{"location":"library.html","page":"Library","title":"Library","text":"Modules = [ManifoldDiff]\nPages = [\"differentials.jl\"]\nOrder = [:type, :function, :constant]\nPrivate = true","category":"page"},{"location":"library.html#ManifoldDiff.differential_exp_argument-Tuple{AbstractManifold, Any, Any, Any}","page":"Library","title":"ManifoldDiff.differential_exp_argument","text":"Z = differential_exp_argument(M, p, X, Y)\ndifferential_exp_argument!(M, Z, p, X, Y)\n\ncomputes D_Xexp_pXY (in place of Z). Note that X   T_X(T_pmathcal M) = T_pmathcal M is still a tangent vector.\n\nSee also\n\ndifferential_exp_basepoint, jacobi_field\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.differential_exp_argument_lie_approx-Tuple{AbstractManifold, Any, Any, Any}","page":"Library","title":"ManifoldDiff.differential_exp_argument_lie_approx","text":"differential_exp_argument_lie_approx(M::AbstractManifold, p, X, Y; n)\n\nApproximate differential of exponential map based on Lie group exponential. The formula reads (see Theorem 1.7 of [Helgason1978])\n\nD_X exp_p(X)Y = (mathrmdL_exp_e(X))_eleft(sum_k=0^nfrac(-1)^k(k+1)(operatornamead_X)^k(Y)right)\n\nwhere (operatornamead_X)^k(Y) is defined recursively as (operatornamead_X)^0(Y) = Y, operatornamead_X^k+1(Y) = X operatornamead_X^k(Y).\n\n[Helgason1978]: S. Helgason, Differential Geometry, Lie Groups, and Symmetric Spaces, First Edition. Academic Press, 1978.\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.differential_exp_basepoint-Tuple{AbstractManifold, Any, Any, Any}","page":"Library","title":"ManifoldDiff.differential_exp_basepoint","text":"Z = differential_exp_basepoint(M, p, X, Y)\ndifferential_exp_basepoint!(M, Z, p, X, Y)\n\nCompute D_pexp_p XY (in place of Z).\n\nSee also\n\ndifferential_exp_argument, jacobi_field\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.differential_inverse_retract_argument_fd_approx-Tuple{AbstractManifold, Any, Any, Any}","page":"Library","title":"ManifoldDiff.differential_inverse_retract_argument_fd_approx","text":"differential_inverse_retract_argument_fd_approx(\n    M::AbstractManifold,\n    p,\n    q,\n    X;\n    retr::AbstractRetractionMethod = default_retraction_method(M),\n    invretr::AbstractInverseRetractionMethod = default_inverse_retraction_method(M),\n    h::Real=sqrt(eps(eltype(X))),\n)\n\nApproximate the differential of the inverse retraction invretr using a finite difference formula (see Eq. (16) in [Zimmermann2019]):\n\nfracoperatornameretr^-1_q(operatornameretr_p(hX)) - operatornameretr^-1_q(operatornameretr_p(-hX))2h\n\nwhere h is the finite difference step h, operatornameretr^-1 is the inverse retraction invretr and operatornameretr is the retraction retr.\n\n[Zimmermann2019]: R. Zimmermann, “Hermite interpolation and data processing errors on Riemannian matrix manifolds,” arXiv:1908.05875 [cs, math], Sep. 2019, Available: http://arxiv.org/abs/1908.05875\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.differential_log_argument-Tuple{AbstractManifold, Any, Any, Any}","page":"Library","title":"ManifoldDiff.differential_log_argument","text":"Y = differential_log_argument(M, p, q, X)\ndifferential_log_argument!(M, Y, p, q, X)\n\ncomputes D_qlog_pqX (in place of Y).\n\nSee also\n\ndifferential_log_basepoint, jacobi_field\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.differential_log_basepoint-Tuple{AbstractManifold, Any, Any, Any}","page":"Library","title":"ManifoldDiff.differential_log_basepoint","text":"Y = differential_log_basepoint(M, p, q, X)\ndifferential_log_basepoint!(M, Y, p, q, X)\n\ncomputes D_plog_pqX (in place of Y).\n\nSee also\n\ndifferential_log_argument, jacobi_field\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.differential_shortest_geodesic_endpoint-Tuple{AbstractManifold, Any, Any, Any, Any}","page":"Library","title":"ManifoldDiff.differential_shortest_geodesic_endpoint","text":"Y = differential_shortest_geodesic_endpoint(M, p, q, t, X)\ndifferential_shortest_geodesic_endpoint!(M, Y, p, q, t, X)\n\nCompute D_qγ(tpq)X (in place of Y).\n\nSee also\n\ndifferential_shortest_geodesic_startpoint, jacobi_field\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.differential_shortest_geodesic_startpoint-Tuple{AbstractManifold, Any, Any, Any, Any}","page":"Library","title":"ManifoldDiff.differential_shortest_geodesic_startpoint","text":"Y = differential_shortest_geodesic_startpoint(M, p, q, t, X)\ndifferential_shortest_geodesic_startpoint!(M, Y, p, q, t, X)\n\nCompute D_p γ(tpq)η (in place of Y).\n\nSee also\n\ndifferential_shortest_geodesic_endpoint, jacobi_field\n\n\n\n\n\n","category":"method"},{"location":"library.html","page":"Library","title":"Library","text":"Modules = [ManifoldDiff]\nPages = [\"diagonalizing_projectors.jl\"]\nOrder = [:type, :function, :constant]\nPrivate = true","category":"page"},{"location":"library.html#ManifoldDiff.AbstractProjector","page":"Library","title":"ManifoldDiff.AbstractProjector","text":"abstract type AbstractProjector end\n\nAn abstract type for projectors on a tangent space T_pM for fixed values of p and M. Calling a projector on a tangent vector returns a new tangent vector:\n\n(Π::AbstractProjector)(X) -> Y\n\nProjectors assume that X is a valid vector from T_pM.\n\n\n\n\n\n","category":"type"},{"location":"library.html#ManifoldDiff.CoprojectorOntoVector","page":"Library","title":"ManifoldDiff.CoprojectorOntoVector","text":"CoprojectorOntoVector{TM<:AbstractManifold,TP,TX}\n\nA structure that represents projector onto the subspace of the tangent space at p from manifold M othogonal to vector X of unit norm.\n\nConstructor\n\nCoprojectorOntoVector(M::AbstractManifold, p, X)\n\n\n\n\n\n","category":"type"},{"location":"library.html#ManifoldDiff.ProjectorOntoVector","page":"Library","title":"ManifoldDiff.ProjectorOntoVector","text":"ProjectorOntoVector{TM<:AbstractManifold,TP,TX}\n\nA structure that represents projector onto the subspace of the tangent space at p from manifold M spanned by tangent vector X of unit norm.\n\nConstructor\n\nProjectorOntoVector(M::AbstractManifold, p, X)\n\n\n\n\n\n","category":"type"},{"location":"library.html#ManifoldDiff.diagonalizing_projectors-Tuple{AbstractManifold, Any, Any}","page":"Library","title":"ManifoldDiff.diagonalizing_projectors","text":"diagonalizing_projectors(M::AbstractManifold, p, X)\n\nCompute eigenvalues of the Jacobi operator Y  R(YX)X, where R is the curvature endomorphism, together with projectors onto eigenspaces of the operator. Projectors are objects of subtypes of AbstractProjector.\n\nBy default constructs projectors using the DiagonalizingOrthonormalBasis.\n\n\n\n\n\n","category":"method"},{"location":"library.html#Gradients","page":"Library","title":"Gradients","text":"","category":"section"},{"location":"library.html","page":"Library","title":"Library","text":"Modules = [ManifoldDiff]\nPages = [\"gradients.jl\"]\nOrder = [:type, :function, :constant]\nPrivate = true","category":"page"},{"location":"library.html#ManifoldDiff.grad_distance","page":"Library","title":"ManifoldDiff.grad_distance","text":"grad_distance(M, q, p[, c=2])\ngrad_distance!(M, X, q, p[, c=2])\n\ncompute the (sub)gradient of the distance (default: squared), in place of X.\n\nf(p) = frac1c d^c_mathcal M(p q)\n\nto a fixed point q on the manifold M and c is an integer. The (sub-)gradient reads\n\noperatornamegradf(p) = -d_mathcal M^c-2(p q)log_pq\n\nfor cneq 1 or pneq  q. Note that for the remaining case c=1, p=q, the function is not differentiable. In this case, the function returns the corresponding zero tangent vector, since this is an element of the subdifferential.\n\nOptional\n\nc – (2) the exponent of the distance,  i.e. the default is the squared distance\n\n\n\n\n\n","category":"function"},{"location":"library.html#ManifoldDiff.subgrad_distance","page":"Library","title":"ManifoldDiff.subgrad_distance","text":"subgrad_distance(M, q, p[, c = 1; atol = 0])\nsubgrad_distance!(M, X, q, p[, c = 1; atol = 0])\n\ncompute the subgradient of the distance (in place of X)\n\nf(p) = frac1c d^c_mathcal M(p q)\n\nto a fixed point q on the manifold M and c is an integer. The subgradient reads\n\npartial f(p) = -d_mathcal M^c-2(p q)log_pq\n\nfor cneq 1 or pneq  q. Note that for the remaining case c=1, p=q, the function is not differentiable. In this case, the subgradient is given by a tangent vector at p with norm less than or equal to one.\n\nOptional\n\nc – (1) the exponent of the distance,  i.e. the default is the distance\natol – (0) the tolerance to use when evaluating the distance between p and q.\n\n\n\n\n\n","category":"function"},{"location":"library.html#Jacobi-fields","page":"Library","title":"Jacobi fields","text":"","category":"section"},{"location":"library.html","page":"Library","title":"Library","text":"Modules = [ManifoldDiff]\nPages = [\"Jacobi_fields.jl\"]\nOrder = [:type, :function, :constant]","category":"page"},{"location":"library.html#ManifoldDiff.adjoint_Jacobi_field-Union{Tuple{Tβ}, Tuple{AbstractManifold, Any, Any, Any, Any, Tβ}} where Tβ","page":"Library","title":"ManifoldDiff.adjoint_Jacobi_field","text":"Y = adjoint_Jacobi_field(M, p, q, t, X, β)\nadjoint_Jacobi_field!(M, Y, p, q, t, X, β)\n\nCompute the AdjointJacobiField J along the geodesic γ_pq on the manifold mathcal M with initial conditions (depending on the application) X  T_γ_pq(t)mathcal M and weights β. The result is a vector Y  T_pmathcal M. The main difference to jacobi_field is the, that the input X and the output Y switched tangent spaces. The computation can be done in place of Y.\n\nFor details see jacobi_field\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.jacobi_field-Union{Tuple{Tβ}, Tuple{AbstractManifold, Any, Any, Any, Any, Tβ}} where Tβ","page":"Library","title":"ManifoldDiff.jacobi_field","text":"Y = jacobi_field(M, p, q, t, X, β)\njacobi_field!(M, Y, p, q, t, X, β)\n\ncompute the Jacobi field J along the geodesic γ_pq on the manifold mathcal M with initial conditions (depending on the application) X  T_pmathcal M and weights β. The result is a tangent vector Y from T_γ_pq(t)mathcal M. The computation can be done in place of Y.\n\nSee also\n\nadjoint_Jacobi_field\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.βdifferential_exp_argument-Tuple{Any, Number, Any}","page":"Library","title":"ManifoldDiff.βdifferential_exp_argument","text":"βdifferential_exp_argument(κ,t,d)\n\nweights for the jacobi_field corresponding to the differential of the geodesic with respect to its start point D_X exp_p XY. They are\n\nβ(κ) = begincases\nfracsinh(dsqrt-κ)dsqrt-κtext if κ  0\n1  text if  κ = 0\nfracsin(dsqrtκ)dsqrtκtext if κ  0\nendcases\n\nSee also\n\ndifferential_exp_argument, jacobi_field\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.βdifferential_exp_basepoint-Tuple{Any, Number, Any}","page":"Library","title":"ManifoldDiff.βdifferential_exp_basepoint","text":"βdifferential_exp_basepoint(κ,t,d)\n\nweights for the jacobi_field corresponding to the differential of the geodesic with respect to its start point D_p exp_p X Y. They are\n\nβ(κ) = begincases\ncosh(sqrt-κ)text if κ  0\n1  text if  κ = 0\ncos(sqrtκ) text if κ  0\nendcases\n\nSee also\n\ndifferential_exp_basepoint, jacobi_field\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.βdifferential_log_argument-Tuple{Any, Number, Any}","page":"Library","title":"ManifoldDiff.βdifferential_log_argument","text":"βdifferential_log_argument(κ,t,d)\n\nweights for the JacobiField corresponding to the differential of the logarithmic map with respect to its argument D_q log_p qX. They are\n\nβ(κ) = begincases\nfrac dsqrt-κ sinh(dsqrt-κ)text if κ  0\n1  text if  κ = 0\nfrac dsqrtκ sin(dsqrtκ)text if κ  0\nendcases\n\nSee also\n\ndifferential_log_basepoint, jacobi_field\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.βdifferential_log_basepoint-Tuple{Any, Number, Any}","page":"Library","title":"ManifoldDiff.βdifferential_log_basepoint","text":"βdifferential_log_basepoint(κ,t,d)\n\nweights for the jacobi_field corresponding to the differential of the geodesic with respect to its start point D_p log_p qX. They are\n\nβ(κ) = begincases\n-sqrt-κdfraccosh(dsqrt-κ)sinh(dsqrt-κ)text if κ  0\n-1  text if  κ = 0\n-sqrtκdfraccos(dsqrtκ)sin(dsqrtκ)text if κ  0\nendcases\n\nSee also\n\ndifferential_log_argument, differential_log_argument, jacobi_field\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.βdifferential_shortest_geodesic_startpoint-Tuple{Any, Any, Any}","page":"Library","title":"ManifoldDiff.βdifferential_shortest_geodesic_startpoint","text":"βdifferential_shortest_geodesic_startpoint(κ,t,d)\n\nweights for the jacobi_field corresponding to the differential of the geodesic with respect to its start point D_x g(tpq)X. They are\n\nβ(κ) = begincases\nfracsinh(d(1-t)sqrt-κ)sinh(dsqrt-κ)\ntext if κ  0\n1-t  text if  κ = 0\nfracsin((1-t)dsqrtκ)sinh(dsqrtκ)\ntext if κ  0\nendcases\n\nDue to a symmetry argument, these are also used to compute D_q g(t pq)η\n\nSee also\n\ndifferential_shortest_geodesic_endpoint, differential_shortest_geodesic_startpoint, jacobi_field\n\n\n\n\n\n","category":"method"},{"location":"library.html#Riemannian-differentials","page":"Library","title":"Riemannian differentials","text":"","category":"section"},{"location":"library.html","page":"Library","title":"Library","text":"Modules = [ManifoldDiff]\nPages = [\"riemannian_diff.jl\"]\nOrder = [:type, :function, :constant]","category":"page"},{"location":"library.html#ManifoldDiff.AbstractRiemannianDiffBackend","page":"Library","title":"ManifoldDiff.AbstractRiemannianDiffBackend","text":"AbstractRiemannianDiffBackend\n\nAn abstract type for backends for differentiation.\n\n\n\n\n\n","category":"type"},{"location":"library.html#ManifoldDiff.RiemannianProjectionBackend","page":"Library","title":"ManifoldDiff.RiemannianProjectionBackend","text":"RiemannianProjectionBackend <: AbstractRiemannianDiffBackend\n\nThis backend computes the differentiation in the embedding, which is currently limited to the gradient. Let mathcal M denote a manifold embedded in some R^m, where m is usually (much) larger than the manifold dimension. Then we require three tools\n\nA function f ℝ^m  ℝ such that its restriction to the manifold yields the cost function f of interest.\nA project function to project tangent vectors from the embedding (at T_pℝ^m) back onto the tangent space T_pmathcal M. This also includes possible changes of the representation of the tangent vector (e.g. in the Lie algebra or in a different data format).\nA change_representer for non-isometrically embedded manifolds, i.e. where the tangent space T_pmathcal M of the manifold does not inherit the inner product from restriction of the inner product from the tangent space T_pℝ^m of the embedding\n\nsee also riemannian_gradient and [AbsilMahonySepulchre2008], Section 3.6.1 for a derivation on submanifolds.\n\n[AbsilMahonySepulchre2008]: Absil, P.-A., Mahony, R. and Sepulchre R., Optimization Algorithms on Matrix Manifolds Princeton University Press, 2008, doi: 10.1515/9781400830244 open access\n\n\n\n\n\n","category":"type"},{"location":"library.html#ManifoldDiff.TangentDiffBackend","page":"Library","title":"ManifoldDiff.TangentDiffBackend","text":"TangentDiffBackend <: AbstractRiemannianDiffBackend\n\nA backend that uses tangent spaces and bases therein to derive an intrinsic differentiation scheme.\n\nSince it works in tangent spaces at argument and function value, methods might require a retraction and an inverse retraction as well as a basis.\n\nIn the tangent space itself, this backend then employs an (Euclidean) AbstractDiffBackend\n\nConstructor\n\nTangentDiffBackend(diff_backend)\n\nwhere diff_backend is an AbstractDiffBackend to be used on the tangent space.\n\nWith the keyword arguments\n\nretraction an AbstractRetractionMethod (ExponentialRetraction by default)\ninverse_retraction an AbstractInverseRetractionMethod LogarithmicInverseRetraction by default)\nbasis_arg an AbstractBasis (DefaultOrthogonalBasis by default)\nbasis_val an AbstractBasis (DefaultOrthogonalBasis by default)\n\n\n\n\n\n","category":"type"},{"location":"library.html#ManifoldDiff.differential-Tuple{AbstractManifold, Any, Real, ManifoldDiff.AbstractRiemannianDiffBackend}","page":"Library","title":"ManifoldDiff.differential","text":"differential(M::AbstractManifold, f, t::Real, backend::AbstractDiffBackend)\ndifferential!(M::AbstractManifold, f, X, t::Real, backend::AbstractDiffBackend)\n\nCompute the Riemannian differential of a curve f ℝto M on a manifold M represented by function f at time t using the given backend. It is calculated as the tangent vector equal to mathrmdf_t(t)1.\n\nThe mutating variant computes the differential in place of X.\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.gradient-Tuple{AbstractManifold, Any, Any, ManifoldDiff.AbstractRiemannianDiffBackend}","page":"Library","title":"ManifoldDiff.gradient","text":"gradient(M::AbstractManifold, f, p, backend::AbstractRiemannianDiffBackend)\ngradient!(M::AbstractManifold, f, X, p, backend::AbstractRiemannianDiffBackend)\n\nCompute the Riemannian gradient f(p) of a real-valued function fmathcal M to ℝ at point p on the manifold M using the specified AbstractRiemannianDiffBackend.\n\nThe mutating variant computes the gradient in place of X.\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.gradient-Tuple{AbstractManifold, Any, Any, ManifoldDiff.TangentDiffBackend}","page":"Library","title":"ManifoldDiff.gradient","text":"gradient(M, f, p, backend::TangentDiffBackend)\n\nThis method uses the internal backend.diff_backend (Euclidean) on the function\n\n    f(operatornameretr_p(cdot))\n\nwhich is given on the tangent space. In detail, the gradient can be written in terms of the backend.basis_arg. We illustrate it here for an AbstractOrthonormalBasis, since that simplifies notations:\n\noperatornamegradf(p) = operatornamegradf(p) = sum_i=1^d g_p(operatornamegradf(p)X_i)X_i\n\t= sum_i=1^d Df(p)X_iX_i\n\nwhere the last equality is due to the definition of the gradient as the Riesz representer of the differential.\n\nIf the backend is a forward (or backward) finite difference, these coefficients in the sum can be approximates as\n\nDF(p)Y  frac1hbigl( f(exp_p(hY)) - f(p) bigr)\n\nwriting p=exp_p(0) we see that this is a finite difference of fcircexp_p, i.e. for a function on the tangent space, so we can also use other (Euclidean) backends\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.hessian-Tuple{AbstractManifold, Any, Any, ManifoldDiff.TangentDiffBackend}","page":"Library","title":"ManifoldDiff.hessian","text":"hessian(M::AbstractManifold, f, p, backend::TangentDiffBackend)\n\nCompute the Hessian of function f at point p using the given backend. The formula for normal coordinate systems from[SommerFletcherPennec2020] is used.\n\n[SommerFletcherPennec2020]: S. Sommer, T. Fletcher, and X. Pennec, “1 - Introduction to differential and Riemannian geometry,” in Riemannian Geometric Statistics in Medical Image Analysis, X. Pennec, S. Sommer, and T. Fletcher, Eds. Academic Press, 2020, pp. 3–37. doi: 10.1016/B978-0-12-814725-2.00008-X.\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.riemannian_Hessian-Tuple{AbstractManifold, Any, Any, Any, Any}","page":"Library","title":"ManifoldDiff.riemannian_Hessian","text":"riemannian_Hessian(M, p, eG, eH, X)\nriemannian_Hessian!(M, Y, p, eG, eH, X)\n\nConvert the Euclidean Hessian eH=operatornameHess tilde f(p) X of a function f colon mathcal M to mathbb R, which is the restriction of tilde f to mathcal M, given additionally the (Euclidean) gradient operatornamegrad tilde f(p).\n\nBy default it uses the following computation:\n\nLet P_X(Y) = D_p(operatornameproj_T_pmathcal M(Y)X, where Y in T_pmathcal M denote the differential of the projection onto T_pmathcal M with respect to p (at X in direction Y). Here X is a tangent vector in the embedding, cf. differential_project_basepoint(M, p, X, Y). We further denote by operatornameproj_N_pmathcal M(X) = X - operatornameproj_T_pmathcal M(X) the projection onto the normal space at p\n\nThen, following Boumal, 2023, Section 5.11 we can compute\n\noperatornameHess f(p)X\n= operatornameproj_T_pmathcal Mbigl(operatornameHess tilde f(p)X)\n+ mathcal P_XBigl( operatornameproj_N_pmathcal Mbigl( operatornamegrad tilde f (p) bigr) Bigr)\n\nThis method can also be implemented directly, if a more efficient/stable version is known.\n\nThe function is inspired by ehess2rhess in the Matlab package Manopt.\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.riemannian_gradient-Tuple{AbstractManifold, Any, Any}","page":"Library","title":"ManifoldDiff.riemannian_gradient","text":"riemannian_gradient(M, p, Y; embedding_metric=EuclideanMetric())\nriemannian_gradient!(M, X, p, Y; embedding_metric=EuclideanMetric())\n\nFor a given gradient Y = operatornamegrad tilde f(p) in the embedding of a manifold, this function computes the Riemannian gradient operatornamegrad f(p) of the function tilde f restricted to the manifold M. This can also be done in place of X.\n\nBy default it uses the following computation: Let the projection Z = operatornameproj_T_pmathcal M(Y) of Y onto the tangent space at p be given, that is with respect to the inner product in the embedding. Then\n\nZ-Y W = 0 text for all  W in T_pmathcal M\n\nor rearranged YW = ZW. We then use the definition of the Riemannian gradient\n\noperatornamegrad f(p) W_p = Df(p)X = operatornamegradf(p) W = operatornameproj_T_pmathcal M(operatornamegradf(p))W\nquadtextfor all  W in T_pmathcal M\n\nComparing the first and the last term, the remaining computation is the function change_representer.\n\nThis method can also be implemented directly, if a more efficient/stable version is known.\n\nThe function is inspired by egrad2rgrad in the Matlab package Manopt.\n\n\n\n\n\n","category":"method"},{"location":"library.html#Manifold-specific-specializations","page":"Library","title":"Manifold-specific specializations","text":"","category":"section"},{"location":"library.html","page":"Library","title":"Library","text":"Modules = [ManifoldDiff]\nPages = [\"manifolds.jl\"]\nOrder = [:type, :function, :constant]","category":"page"},{"location":"library.html#Differentiation-backends","page":"Library","title":"Differentiation backends","text":"","category":"section"},{"location":"library.html#EmbeddedDiff","page":"Library","title":"EmbeddedDiff","text":"","category":"section"},{"location":"library.html","page":"Library","title":"Library","text":"Modules = [ManifoldDiff]\nPages = [\"embedded_diff.jl\"]\nOrder = [:type, :function, :constant]","category":"page"},{"location":"library.html#ManifoldDiff.ExplicitEmbeddedBackend","page":"Library","title":"ManifoldDiff.ExplicitEmbeddedBackend","text":"ExplicitEmbeddedBackend{TF<:NamedTuple} <: AbstractDiffBackend\n\nA backend to use with the RiemannianProjectionBackend or the TangentDiffBackend, when you have explicit formulae for the gradient in the embedding available.\n\nConstructor\n\nExplicitEmbeddedBackend(M::AbstractManifold; kwargs)\n\nConstruct an ExplicitEmbeddedBackend in the embedding M, where currently the following keywords may be used\n\ngradient for a(n allocating) gradient function gradient(M, p) defined in the embedding\ngradient! for a mutating gradient function gradient!(M, X, p).\n\nNote that the gradient functions are defined on the embedding manifold M passed to the Backend as well\n\n\n\n\n\n","category":"type"},{"location":"library.html#ForwardDiff.jl","page":"Library","title":"ForwardDiff.jl","text":"","category":"section"},{"location":"library.html","page":"Library","title":"Library","text":"Modules = [ManifoldDiff]\nPages = [\"forward_diff.jl\"]\nOrder = [:type, :function, :constant]","category":"page"},{"location":"library.html#ManifoldDiff.ForwardDiffBackend","page":"Library","title":"ManifoldDiff.ForwardDiffBackend","text":"ForwardDiffBackend <: AbstractDiffBackend\n\nDifferentiation backend based on the ForwardDiff.jl package.\n\n\n\n\n\n","category":"type"},{"location":"library.html#FiniteDifferenes.jl","page":"Library","title":"FiniteDifferenes.jl","text":"","category":"section"},{"location":"library.html","page":"Library","title":"Library","text":"Modules = [ManifoldDiff]\nPages = [\"finite_differences.jl\"]\nOrder = [:type, :function, :constant]","category":"page"},{"location":"library.html#ManifoldDiff.FiniteDifferencesBackend","page":"Library","title":"ManifoldDiff.FiniteDifferencesBackend","text":"FiniteDifferencesBackend(method::FiniteDifferenceMethod = central_fdm(5, 1))\n\nDifferentiation backend based on the FiniteDifferences.jl package.\n\n\n\n\n\n","category":"type"},{"location":"library.html#Internal-functions","page":"Library","title":"Internal functions","text":"","category":"section"},{"location":"library.html","page":"Library","title":"Library","text":"Modules = [ManifoldDiff]\nPages = [\"ManifoldDiff.jl\"]\nOrder = [:type, :function, :constant]\nPrivate = true\nPublic=false","category":"page"},{"location":"library.html#ManifoldDiff.AbstractDiffBackend","page":"Library","title":"ManifoldDiff.AbstractDiffBackend","text":"AbstractDiffBackend\n\nAn abstract type for diff backends. See FiniteDifferencesBackend for an example.\n\n\n\n\n\n","category":"type"},{"location":"library.html#ManifoldDiff.CurrentDiffBackend","page":"Library","title":"ManifoldDiff.CurrentDiffBackend","text":"CurrentDiffBackend(backend::AbstractDiffBackend)\n\nA mutable struct for storing the current differentiation backend in a global constant _current_default_differential_backend.\n\nSee also\n\nAbstractDiffBackend, default_differential_backend, set_default_differential_backend!\n\n\n\n\n\n","category":"type"},{"location":"library.html#ManifoldDiff._derivative","page":"Library","title":"ManifoldDiff._derivative","text":"_derivative(f, t[, backend::AbstractDiffBackend])\n\nCompute the derivative of a callable f at time t computed using the given backend, an object of type AbstractDiffBackend. If the backend is not explicitly specified, it is obtained using the function default_differential_backend.\n\nThis function calculates plain Euclidean derivatives, for Riemannian differentiation see for example differential.\n\nnote: Note\nNot specifying the backend explicitly will usually result in a type instability and decreased performance.\n\n\n\n\n\n","category":"function"},{"location":"library.html#ManifoldDiff._gradient","page":"Library","title":"ManifoldDiff._gradient","text":"_gradient(f, p[, backend::AbstractDiffBackend])\n\nCompute the gradient of a callable f at point p computed using the given backend, an object of type AbstractDiffBackend. If the backend is not explicitly specified, it is obtained using the function default_differential_backend.\n\nThis function calculates plain Euclidean gradients, for Riemannian gradient calculation see for example gradient.\n\nnote: Note\nNot specifying the backend explicitly will usually result in a type instability and decreased performance.\n\n\n\n\n\n","category":"function"},{"location":"library.html#ManifoldDiff._hessian","page":"Library","title":"ManifoldDiff._hessian","text":"_hessian(f, p[, backend::AbstractDiffBackend])\n\nCompute the Hessian of a callable f at point p computed using the given backend, an object of type AbstractDiffBackend. If the backend is not explicitly specified, it is obtained using the function default_differential_backend.\n\nThis function calculates plain Euclidean Hessian.\n\nnote: Note\nNot specifying the backend explicitly will usually result in a type instability and decreased performance.\n\n\n\n\n\n","category":"function"},{"location":"library.html#ManifoldDiff._jacobian","page":"Library","title":"ManifoldDiff._jacobian","text":"_jacobian(f, p[, backend::AbstractDiffBackend])\n\nCompute the Jacobian of a callable f at point p computed using the given backend, an object of type AbstractDiffBackend. If the backend is not explicitly specified, it is obtained using the function default_differential_backend.\n\nThis function calculates plain Euclidean Jacobians, for Riemannian Jacobian calculation see for example gradient.\n\nnote: Note\nNot specifying the backend explicitly will usually result in a type instability and decreased performance.\n\n\n\n\n\n","category":"function"},{"location":"library.html#ManifoldDiff.default_differential_backend-Tuple{}","page":"Library","title":"ManifoldDiff.default_differential_backend","text":"default_differential_backend() -> AbstractDiffBackend\n\nGet the default differentiation backend.\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff.set_default_differential_backend!-Tuple{ManifoldDiff.AbstractDiffBackend}","page":"Library","title":"ManifoldDiff.set_default_differential_backend!","text":"set_default_differential_backend!(backend::AbstractDiffBackend)\n\nSet current backend for differentiation to backend.\n\n\n\n\n\n","category":"method"},{"location":"library.html#ManifoldDiff._current_default_differential_backend","page":"Library","title":"ManifoldDiff._current_default_differential_backend","text":"_current_default_differential_backend\n\nThe instance of CurrentDiffBackend that stores the globally default differentiation backend.\n\n\n\n\n\n","category":"constant"},{"location":"index.html#ManifoldDiff","page":"Home","title":"ManifoldDiff","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"The package ManifoldDiff aims to provide automatic calculation of Riemannian gradients of functions defined on manifolds. It builds upon Manifolds.jl.","category":"page"},{"location":"index.html#Naming-scheme","page":"Home","title":"Naming scheme","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Providing a derivative, differential or gradient for a given function, this package adds that information to the function name. For example","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"grad_f for a gradient operatornamegrad f\nsubgrad_f for a subgradient partial f\ndifferential_f for Df (also called pushforward)\ndifferential_f_variable if f has multiple variables / parameters, since a usual writing in math is f_x in this case\nadjoint_differential_f for pullbacks\nadjoint_differential_f_variable if f has multiple variables / parameters\nf_derivative for f","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"the scheme is not completely fixed but tries to follow the mathematical notation.","category":"page"}]
}
