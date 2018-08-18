const OnePoint{N1,T1} = Union{Vector{T1}, NTuple{N1,T1}} where {N1,T1}
const ValueSeq{N2,T2} = Union{Vector{T3}, NTuple{N2,T2}} where {N2,T2}
const PointSeq{N3,T3} = Union{Vector{T3}, NTuple{N3,T3}} where {N3,T2}
