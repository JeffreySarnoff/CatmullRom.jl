root4(x) = sqrt(sqrt(x))

distancesquared(x::T, y::T) where {T<:Number} = x*x + y*y
distancesquared(x::T, y::T, z::T) where {T<:Number} = x*x + y*y + z*z
distancesquared(w::T, x::T, y::T, z::T) where {T<:Number} = w*w + x*x + y*y + z*z
distance(x::T, y::T) where {T<:Number} = sqrt(distancesquared(x,y))
distance(x::T, y::T, z::T) where {T<:Number} = sqrt(distancesquared(x,y,z))
distance(w::T, x::T, y::T, z::T) where {T<:Number} = sqrt(distancesquared(w,x,y,z))

distancesquared(a::A, b::A) where {N,T,S,A<:Union{NTuple{N,T},AbstractVector{S}}} = sum(subsq(a, b))

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
