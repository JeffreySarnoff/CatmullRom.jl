#=
   Determine two extremal points (outside of the boundry points)
   for use in anchoring the sequence of centripetal Catmull-Rom
   spline spans to incorporate the boundry points in the spline.
=#

function outside(initialpoints::T, finalpoints::T) where {T}
    ninitials = length(initialpoints)
    nfinals   = length(finalpoints)
    (ninitials > 1 && nfinals > 1) || throw(DomainError("missing points"))
    
    xinitial = reflectback(initialpoints[1:2]...,)[1]
    xfinal   = reflectforward(finalpoints[end-1:end]...,)[1]
    initial = ninitials > 2 ? thiele3(initialpoints[1:3], xinitial) : linear(initialpoints[1:2], xinitial)
    final  = ninitials > 2 ? thiele3(finalpoints[1:3], xfinal) : linear(finalpoints[1:2], xfinal)

    return initial, final
end
