function Base.push!(collection::NTuple{N,NTuple{M,T}}, item::NTuple{M,T}) where {N,M,T}
    return (collection..., item)
end

function Base.push!(collection::NTuple{N1,NTuple{M,T}}, items::NTuple{N2,NTuple{M,T}}) where {N1,N2,M,T}
    return (collection..., items...)
end

function Base.pushfirst!(collection::NTuple{N,NTuple{M,T}}, item::NTuple{M,T}) where {N,M,T}
    return (item, collection...)
end

function Base.pushfirst!(collection::NTuple{N1,NTuple{M,T}}, items::NTuple{N2,NTuple{M,T}}) where {N1,N2,M,T}
    return (items..., collection...)
end

function myunzip(x)
    nd = length(x[1])
    npts = length(x)
    typ = eltype(eltype(x[1]))
    res = Array{typ, 2}(undef, (npts,nd))
    for k=1:npts
       for i=1:nd
          res[k,i] = x[k][i]
       end
    end
    return res
end

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
