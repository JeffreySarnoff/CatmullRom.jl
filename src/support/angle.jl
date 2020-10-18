
#=
  relatively fast determination of angular separation
     UNCHECKED PRECONDITION:
        both points are given relative to the same origin
 
   >>>  for a numerically rigourous approach to angular separation
   >>>  use AngleBetweenVectors.jl
=#

function anglesep(pointa::P1, pointb::P1) where {P1}
    dota = dot(pointa, pointa)
    dotb = dot(pointb, pointb)
    (iszero(dota) || iszero(dotb)) && return zero(T)

    dotb = sqrt(dota * dotb)
    dota = dot(pointa, pointb)
    acos( dota / dotb )
end
