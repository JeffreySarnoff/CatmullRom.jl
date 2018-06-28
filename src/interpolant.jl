
"""
    into01((xs...,))
    into01([xs...,])

maps values into 0.0:1.0 (minimum(xs) --> 0.0, maximum(xs) --> 1.0)
"""
function into01(values::U) where {N,T, U<:Union{NTuple{N,T}, Vector{T}}}
    mn, mx = minimum(values), maximum(values)
    delta = mx - mn
    map(clamp01, (values .- mn) ./ delta)
end

@inline clamp01(x::T) where {T<:Real} = clamp(x, zero(T), one(T))
