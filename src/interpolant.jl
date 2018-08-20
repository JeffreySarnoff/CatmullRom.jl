"""
    clamp01(x)

min(max(0, x), 1)
"""
@inline clamp01(x::T) where {T<:Real} = clamp(x, zero(T), one(T))

"""
    within01(n::Int)

generate n >= 0 values uniformly spaced inside 0..1, excludes 0, 1

# Examples
    within01(0) == []
    within01(1) == [1/2]
    within01(2) == [1/3, 2/3]
"""
function within01(n::Int)
    n = flipsign(n,n)                   # abs(n)
    nplus1 = n + 1
    inv_nplus1 = inv(nplus1)
    collect(inv_nplus1:inv_nplus1:n/nplus1)
end

"""
    uniform01(n::Int)

generate n >= 2 values uniformly spaced from 0..1, includes 0, 1

# Examples
    uniform01(2) == [0.0, 1.0]
    uniform01(3) == [0.0, 0.5, 1.0]
    uniform01(4) == [0.0, 1/3, 2/3, 1.0]
"""
uniform01(n::Int) = collect(0.0:inv(n+1):1.0)


#!!REVIEW IMPLEMENTATION!!

"""
    into01((xs...,))

maps values into 0.0:1.0 (minimum(xs) --> 0.0, maximum(xs) --> 1.0)
"""
function into01(unitrange::UnitRange{T}) where {T}
    span = unitrange.stop - unitrange.start
    return collect(LinRange(0.0, 1.0, span))
end

function into01(values::U) where {U<:Union{AbstractVector{T}, NTuple{N,T}}} where {N,T}
    valuemin = minimum(values)
    span = maximum(values) - valuemin
    delta = sqrt(sum(map(x->x*x, span)))
    result = collect(clamp01((v .- valuemin)./span) for v in values)
    return result
end


