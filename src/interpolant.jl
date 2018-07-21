"""
    into01((xs...,))
    into01([xs...,])

maps values into 0.0:1.0 (minimum(xs) --> 0.0, maximum(xs) --> 1.0)
"""
function into01(values::U) where {N,T, U<:Union{NTuple{N,T}, Vector{T}}}
    mn, mx = minimum(values), maximum(values)
    delta = mx .- mn
    delta = sqrt(sum(map(x->x*x, delta)))
    result = collect((v .- mn)./delta for v in values)
    #result = map(clamp01, result)
    return result
end

@inline clamp01(x::T) where {T<:Real} = clamp(x, zero(T), one(T))

# uniform separation in 0..1 inclusive
uniformsep(n::Int) = n >= 2 ? collect(0.0:inv(n-1):1.0) : throw(DomainError("$n < 2"))

# roots of the Chebyshev polynomials of the second kind, mapped into 0..1

chebroot2_01(k::Int, n::Int) = (1-cospi((k-1)/(n-1)))/2
chebroots2_01(n::Int) = n >= 1 ? [chebroot2_01(k,n) for k=1:n] : throw(DomainError("$n < 1"))

# Chebyshev type 2 roots mapped into 0..1 with 0 and 1 included
"""
    chebyshevsep(n)


Chebyshev type 2 roots mapped into 0..1, with 0 and 1 included

chebyshevsep(5) gives five values: [0.0, b, c, d, 1.0]
"""
chebyshevsep(n::Int) = n >= 0 ? chebroots2_01(n) : throw(DomainError("$n < 0"))
