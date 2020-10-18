#= 
`AbstractPoint{ND,T}` is an abstraction over all internally recognizable point representations.

With any of the available parameter-realized sorts of _Point_ (`Point{ND,T}`) 
- the parameter `ND` is to be assigned a strictly postive Int that reflects the dimensionality
- the parameter `T`  is to be assigned the numeric type that carries these coordinate values
=#

abstract type AbstractPoint{ND, T} end

struct NTupPoint{ND, T} <: AbstractPoint{ND, T}
    coords::NTuple{ND, T}
end

struct SVecPoint{ND, T} <: AbstractPoint{ND, T}
    coords::SVector{ND, T}
end

const Point{ND,T} = Union{NTupPoint{ND,T}, SVecPoint{ND,T}}

Point(coords::NTuple{ND,T}) where {T,ND} = NTupPoint{ND,T}(coords)
Point(coords::SVector{ND,T}) where {T,ND} = SVecPoint{ND,T}(coords)
Point(coords::Vector{T}) where {T} = SVecPoint{length(coords),T}(coords)

Base.convert(::Type{SVecPoint{ND,T}, x::NTupPoint{ND,T}) where {T,ND} = SVecPoint{ND,T}(x.coords)
Base.convert(::Type{NTupPoint{ND,T}, x::SVecPoint{ND,T}) where {T,ND} = NTupPoint{ND,T}(x.coords.data)
    
Base.promote_rule(::Type{NTupPoint{ND,T}}, ::Type{SVecPoint{ND,T}}) where {T,ND} = SVecPoint{ND,T}
