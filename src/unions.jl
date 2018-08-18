const OnePoint{N,T} = NTuple{N,T} where {N,T}
const Items{N,T} = Union{A, NT} where {N,T, A<:Vector{T}, NT<:NTuple{N,T}}
