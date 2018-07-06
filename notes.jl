julia> cmrpoly(u) = [1 u u^2 u^3]
cmrpoly (generic function with 2 methods)

julia> cmrpoints(p1,p2,p3,p4) = Matrix(transpose([p1 p2 p3 p4]))
cmrpoints (generic function with 1 method)

julia> cmrmatrix(τ) = Matrix(transpose(reshape([0 1 0 0 -τ 0 τ 0 2τ τ-3 3-2τ -τ -τ 2-τ τ-2 τ], (4,4))))
cmrmatrix (generic function with 1 method)

julia> cmrpoly(1//7) * cmrmatrix(1//2) * cmrpoints([0,1],[1,5],[2,4],[3,2])
1×2 Array{Rational{Int64},2}:
 8//7  1759//343

