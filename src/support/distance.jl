root4(x) = sqrt(sqrt(x))

distancesquared(a::T, b::T) where {T} = sum(subsq(a, b))

subsq(a::NTuple{N,T}, b::NTuple{N,T}) where {N,T} =
    (((a[i] - b[i])^2 for i=1:N)...,)

subsq(a::T, b::T) where {S,T<:AbstractVector{S}} =
    (((a[i] - b[i])^2 for i=1:min(length(a),length(b)))...,)

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
