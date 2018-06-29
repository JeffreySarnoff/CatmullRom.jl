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

# uniform separation in 0..1 inclusive
uniformsep(n::Int) = n >= 2 ? collect(0.0:inv(n-1):1.0) : throw(DomainError("$n < 2"))

# Chebyshev type 1 roots mapped into 0..1
chebroot01(k::Int, n::Int) = (1 + cospi( (2*(n-k)+1) / (2*n) )) / 2
chebroots01(n::Int) = n >= 1 ? [chebroot01(k,n) for k=1:n] : throw(DomainError("$n < 1"))

# Chebyshev type 1 roots mapped into 0..1 with 0 and 1 appended
"""
    zero_chebroots_one(n)


Chebyshev type 1 roots mapped into 0..1, with 0 and 1 appended
"""
zero_chebroots_one(n::Int) = n >= 0 ? [0.0, chebroots01(n)..., 1.0] : throw(DomainError("$n < 0"))


function validate_interpolants(interpolants::Union{A,NTuple{N,F}}) where {A<:AbstractArray, N, D, T, F}
   interpolants[1] == 0 && interpolants[end] == 1 &&
       minimum(interpolants) == 0 && maximum(interpolants) == 1 && 
       length(interpolants) == length(unique(interpolants)) ||
       throw(ErrorException("improper interpolants: $interpolants"))
end

function fixup(interpolants::U2) where (U2)
    interps = sort([interpolants...,])
    if interps[end] < 1.0
        interps = [interps..., 1.0]
    end
    if 0.0 < interps[1] < 1.0
        interps = [0.0, interps...,]
    end
    if interps[1] != 0 || interps[end] != 1
        interps = into01(interps)
    end
    return (interps...,)
end

# mapping into centripetal curve 

sqr(x) = x * x
distance(pa, pb) = sqrt(sum(sqr.(pb .- pa)))
dist(pa, pb) = sqrt(distance(pa, pb))

# dists = [dist(pts[i],pts[i+1]) for i=1:length(pts)-1]

dists(pts) = [dist(pts[i,:],pts[i+1,:]) for i=1:length(pts[:,1])-1]

#totaldist = sum(dists)
centripetals(ptdists) = [0.0, (cumsum(ptdists) ./ sum(ptdists))...,]

centripetals(pts) = centripetals(dists(pts))


