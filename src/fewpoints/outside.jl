function outside(points)
    initial, final = outside(points[1:4], points[end-3:end])
    return [initial, points..., final]
end

#=
   Determine two extremal points (outside of the boundry points)
   for use in anchoring the sequence of centripetal Catmull-Rom
   spline spans to incorporate the boundry points in the spline.
=#

function outside(initialpoints::T, finalpoints::T) where {T}
    ninitials = length(initialpoints)
    nfinals   = length(finalpoints)
    (ninitials > 1 && nfinals > 1) || throw(DomainError("missing points"))
    
    initialpt = reflectback(initialpoints[1:2]...,)
    finalpt   = reflectforward(finalpoints[end-1:end]...,)
    initial = ninitials > 3 ? thiele4(initialpoints[1:4], initialpt[1]) : 
              ninitials > 2 ? thiele3(initialpoints[1:3], initialpt[1]) : initialpt
    final   = nfinals > 3 ? thiele4(finalpoints[1:4], finalpt[1]) : 
              nfinals > 2 ? thiele3(finalpoints[1:3], finalpt[1]) : finalpt
    return initial, final
end
