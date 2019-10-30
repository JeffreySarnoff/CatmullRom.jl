vector(x::Vector{T}) where {T} = x
vector(x::Array{T,N}) where {T,N} = vec(x)

function vector(x::NTuple{N,T}) where {N,T}
     result = Vector{T}(undef, N)
     for i=1:N
         result[i] = x[i]
     end
     return result
end
    

"""
     coordvecs(seq_of_coord, ::Type{T}=Float64)

``` 
vec_of_coord_tuples = [(1, 2), (2, 4)]
coordvecs(vec_of_coord_tuples) == [[1.0, 2.0], [2.0, 4.0]]

vec_of_coord_vectors = [[1, 2], [2, 4]]
coordvecs(vec_of_coord_vectors) == [[1.0, 2.0], [2.0, 4.0]]

tuple_of_coord_tuples = ((1, 2), (2, 4))
coordvecs(tuple_of_coord_tuples) == ([1.0, 2.0], [2.0, 4.0])
coordvecs(tuple_of_coord_tuples, Float32) == ([1.0f0, 2.0f0], [2.0f0, 4.0f0])
```
"""
coordvecs(coordseq::Vector{Vector{T}}) where {T<:Integer} = coordvecs(coordseq, Float64)
coordvecs(coordseq::Vector{NTuple{N,T}}) where {N, T<:Integer} = coordvecs(coordseq, Float64)
coordvecs(coordseq::NTuple{N1, NTuple{N2,T}}) where {N1, N2, T<:Integer} = coordvecs(coordseq, Float64)
coordvecs(coordseq::NTuple{N, Vector{T}}) where {N, T<:Integer} = coordvecs(coordseq, Float64)

coordvecs(coordseq::Vector{Vector{T}}) where {T<:AbstractFloat} = coordvecs(coordseq, T)
coordvecs(coordseq::Vector{NTuple{N,T}}) where {N, T<:AbstractFloat} = coordvecs(coordseq, T)
coordvecs(coordseq::NTuple{N1, NTuple{N2,T}}) where {N1, N2, T<:AbstractFloat} = coordvecs(coordseq, T)
coordvecs(coordseq::NTuple{N, Vector{T}}) where {N, T<:AbstractFloat} = coordvecs(coordseq, T)

coordvecs(coordseq::Vector{Vector{T1}}, ::Type{T2}) where {T1<:Real, T2<:AbstractFloat} = map(x->T2.(x), coordseq)
coordvecs(coordseq::Vector{NTuple{N,T1}}, ::Type{T2}) where {N, T1<:Real, T2<:AbstractFloat} = map(x->[T2.(x)...,], coordseq)
coordvecs(coordseq::NTuple{N1, NTuple{N2,T1}}, ::Type{T2}) where {N1, N2, T1<:Real, T2<:AbstractFloat} = map(x->[T2.(x)...,], coordseq)
coordvecs(coordseq::NTuple{N, Vector{T1}}, ::Type{T2}) where {N, T1<:Real, T2<:AbstractFloat} = map(x->[T2.(x)...,], coordseq)

vectorofcoordvecs(coordseq) = vector(coordvecs(coordseq))
