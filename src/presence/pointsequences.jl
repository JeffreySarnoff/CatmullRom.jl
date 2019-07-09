# The sorts of structures understood to hold the coordinates of a point
const NumVec = Array{T,1} where {T<:Number};
const NumTup = NTuple{N,T} where {N,T<:Number};
const NumNT = NamedTuple{S,NTuple{N,T}} where {S,N,T<:Number};
const OnePoint = Union{NumVec, NumTup, NumNT}

# The sorts of sequences understood to hold point coordinates
const VecNumVec = AbstractArray{Array{T,1},1} where {T<:Number};
const VecNumTup = AbstractArray{NTuple{N,T},1} where {N,T<:Number};
const TupNumTup = NTuple{M,NTuple{N,T}} where {M,N,T<:Number};
const TupNumVec = NTuple{M,Array{T,1}} where {M,N,T<:Number};
const VecNumNT = AbstractArray{NamedTuple{S,NTuple{N,T}},1} where {S,N,T<:Number};
const TupNumNT = NTuple{M,NamedTuple{S,NTuple{N,T}}} where {M,S,N,T<:Number};

const Points = Union{VecNumVec, VecNumTup, TupNumTup, TupNumVec, VecNumNT, TupNumNT};

isclosed(firstpoint::OnePoint, lastpoint::OnePoint) = firstpoint == lastpoint
isclosed(points::Points) = isclosed(first(points), last(points))

npoints(pts::Points) = length(pts)
npoints(pts::Base.Iterators.Zip) = length(pts)

ncoords(pts::Points) = eltype(pts) <: NamedTuple ? length(Tuple(pts[1])) : length(pts[1])
ncoords(pts::Base.Iterators.Zip) = eltype(pts) <: NamedTuple ? length(Tuple(first(pts))) : length(first(pts))

coordtype(pts::Points) = eltype(eltype(pts))
coordtype(pt::OnePoint) = eltype(pt)

coordtype(x::T) where {T} = coordtype(T)

function coordtype(::Type{T}) where {T}
    result = T
    !(result<:Number) ? (result = eltype(result)) : (return result)
    !(result<:Number) ? (result = eltype(result)) : (return result)
    !(result<:Number) ? (result = eltype(result)) : (return result)
    !(result<:Number) ? throw(ErrorException("unable to discern the coordinate type for $T")) : (return result)
end

Base.convert(::Type{Array{T,1}}, x::NTuple{N,T}) where {N,T} = [x...,]
Base.convert(::Type{NTuple{N,T}}, x::Array{T,1}) where {N,T} = (x...,)

function coords_to_cols(coords...)
    n_cols = length(coords)
    n_rows = length(coords[1])
    return reshape(vcat(coords...,), n_rows, n_cols)
end

cols_to_coords(cols::Array{T,2}) where T = [cols[:,i] for i=1:(size(cols)[2])]
