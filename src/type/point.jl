#= 
`AbstractPoint{ND,T}` is an abstraction over all internally recognizable point representions.

With any of the available parameter-realized sorts of _Point_ (`Point{ND,T}`) 
- the parameter `ND` is to be assigned a strictly postive Int that reflects the dimensionality
- the parameter `T`  is to be assigned the numeric type that carries these coordinate values
=#

abstract type AbstractPoint{ND, T} end

struct Point{ND, T} <: AbstractPoint{ND, T}
    coords::NTuple{ND, T}
end
