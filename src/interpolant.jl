@inline clamp01(x::T) where {T<:Real} = clamp(x, zero(T), one(T))

"""
    into01((xs...,))

maps values into 0.0:1.0 (minimum(xs) --> 0.0, maximum(xs) --> 1.0)
"""
function into01(values::U) where {N,T, U<:Union{NTuple{N,T}, Vector{T}}}
    mn, mx = minimum(values), maximum(values)
    delta = mx .- mn
    delta = sqrt(sum(map(x->x*x, delta)))
    result = collect(clamp01((v .- mn)./delta) for v in values)
    return result
end

function prepoint(point1, point2)
    point1 - (point2 - point1)/8
end
function postpoint(point1, point2)
    point2 + (point2 - point1)/8
end

function thiele3(point1, point2, point3, x)
    t1 = point1[2] - point2[2]
    t2 = point2[2] - point3[2]
    t1 = inv(t1)
    t2 = inv(t2)
    t1 = (point1[1] - point2[1]) * t1
    t2 = -(point2[1] - point3[1]) * t2 + t1
    t2 = inv(t2)
    t2 = (point1[1] - point3[1]) * t2 - point1[2] + point2[2]
    t2 = inv(t2_
    t1 = -(point2[1] - x) * t2 + t1
    t1 = inv(t1)
    return (-(point1[1] - x) * t1 + point1[2])
end

# uniform separation in 0..1 inclusive
uniformsep(n::Int) = n >= 2 ? collect(0.0:inv(n-1):1.0) : throw(DomainError("$n < 2"))

# chebyshev T(x) zeros in 0..1 inclusive
chebyshevsep(n::Int) = n >= 2 ? chebTzeros01(n) : throw(DomainError("$n < 2"))

# roots of Chebyshev polynomials (T,U,V,W)
#    shifted from [-1,+1] into [0,+1]

# roots of T(x)

chebTzero(n, k) = chebTzero(Float64, n, k)
chebTzero(::Type{T}, n, k) where {T} = cospi(T(n-k + 1/2) / n)

shift_chebTzero(n, k) = shift_chebTzero(Float64, n, k)
shift_chebTzero(::Type{T}, n, k) where {T} = (chebTzero(n,k) + 1) / 2

chebTzeros(n, k) = chebTzeros(Float64, n, k)
function chebTzeros(::Type{T}, n) where {T}
    result = Vector{T}(undef, n)
    for k = 1:n
        result[k] = chebTzero(T,n,k)
    end
    return result
end

chebTzeros01(n, k) = chebTzeros01(Float64, n, k)
function chebTzeros01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = shift_chebTzero(T,n,k)
    end
    return result
end

# roots of U(x)

chebUzero(n, k) = chebUzero(Float64, n, k)
chebUzero(::Type{T}, n, k) = cospi(T(n-k+1) / T(n+1))

shift_chebUzero(n, k) = shift_chebUzero(Float64, n, k)
shift_chebUzero(::Type{T}, n, k) = (chebUzero(T,n,k) + 1) / 2

chebUzeros(n, k) = chebUzeros(Float64, n, k)
function chebUzeros(::Type{T}, n) where {T}
    result = Vector{T}(undef, n)
    for k = 1:n
        result[k] = shift_chebUzero(T,n,k)
    end
    return result
end

chebUzeros01(n, k) = chebUzeros01(Float64, n, k)
function chebUzeros01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = shift_chebUzero(T,n,k)
    end
    return result
end

# zeros of V(x)

chebVzero(n, k) = chebVzero(Float64, n, k)
chebVzero(::Type{T}, n, k) = cospi(T(n-k+1/2) / T(n+1/2))

shift_chebVzero(n, k) = shift_chebVzero(Float64, n, k)
shift_chebVzero(::Type{T}, n, k) = (chebVzero(T,n,k) + 1) / 2

chebVzeros(n, k) = chebVzeros01(Float64, n, k)
function chebVzeros(::Type{T}, n) where {T}
    result = Vector{T}(undef, n)
    for k = 1:n
        result[k] = shift_chebVzero(T,n,k)
    end
    return result
end

chebVzeros01(n, k) = chebVzeros01(Float64, n, k)
function chebVzeros01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = shift_chebVzero(T,n,k)
    end
    return result
end

# zeros of W(x)


chebWzero(::Type{T}, n, k) = cospi(T(n-k+1) / T(n+1/2))

shift_chebWzero(n, k) = shift_chebWzero(Float64, n, k)
shift_chebWzero(::Type{T}, n, k) = (chebWzero(T,n,k) + 1) / 2

chebWzeros(n, k) = chebWzeros(Float64, n, k)
function chebWzeros(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    for k = 1:n
        result[k] = shift_chebWzero(T,n,k)
    end
    return result
end

chebWzeros01(n, k) = chebWzeros01(Float64, n, k)
function chebWzeros01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = shift_chebWzero(T,n,k)
    end
    return result
end

# extrema of _weighted_ Chebyshev polynomials (types T,U,V,W or 1,2,3,4)
# includes reverse order 1-extrema with extrema 

# extrema of T(x)

chebTextrema(n, k) = chebTextrema(Float64, n, k)
chebTextrema(::Type{T}, n, k) where {T} = cospi(T(n-k) / n)

shift_chebTextrema(n, k) = shift_chebTextrema(Float64, n, k) 
shift_chebTextrema(::Type{T}, n, k) where {T} = (chebTextrema(T, n, k) + 1)/2

chebTextrema(n, k) = chebTextrema(Float64, n, k)
function chebTextrema(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    for k = 1:n
        result[k] = shift_chebTextrema(T,n,k)
    end
    return result
end

chebTextrema01(n, k) = chebTextrema01(Float64, n, k)
function chebTextrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = shift_chebTextrema(T,n,k)
    end
    return result
end



# extrema of sqrt(1-x^2) * U(x)

chebUextrema(n, k) = chebUextrema(Float64, n, k) 
chebUextrema(::Type{T}, n, k) where {T} = cospi(T(2*(n-k)+1) / T(2*(n+1)))

shift_chebUextrema(n, k) = shift_chebUextrema(Float64, n, k) 
shift_chebUextrema(::Type{T}, n, k) where {T} = (chebUextrema(T, n, k) + 1)/2

chebUextrema(n, k) = chebUextrema(Float64, n, k)
function chebUextrema(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    for k = 1:n
        result[k] = shift_chebUextrema(T,n,k)
    end
    return result
end

chebUextrema01(n, k) = chebUextrema01(Float64, n, k)
function chebUextrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = shift_chebUextrema(T,n,k)
    end
    return result
end

chebUUextrema01(n, k) = chebUUextrema01(Float64, n, k)
function chebUUextrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, 2*n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    result[2:n+1] = chebUextrema(T,  n, k)
    result[n+2:(2*n+1)] = one(T) .- result[2:n+1]
    sort!(result)
    return result
end


# extrema of sqrt(1+x) * V(x)

chebVextrema(n, k) = chebVextrema(Float64, n, k) 
chebVextrema(::Type{T}, n, k) where {T} = cospi(T(2*(n-k)) / T(2*n+1))

shift_chebVextrema(n, k) = shift_chebVextrema(Float64, n, k) 
shift_chebVextrema(::Type{T}, n, k) where {T} = (chebVextrema(T, n, k) + 1)/2

chebVextrema(n, k) = chebVextrema(Float64, n, k)
function chebVextrema(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    for k = 1:n
        result[k] = shift_chebVextrema(T,n,k)
    end
    return result
end

chebVextrema01(n, k) = chebVextrema01(Float64, n, k)
function chebVextrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = shift_chebVextrema(T,n,k)
    end
    return result
end

chebVVextrema01(n, k) = chebVVextrema01(Float64, n, k)
function chebVVextrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, 2*n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    result[2:n+1] = chebVextrema(T,  n, k)
    result[n+2:(2*n+1)] = one(T) .- result[2:n+1]
    sort!(result)
    return result
end

# extrema of sqrt(1-x) * W(x)

chebWextrema(n, k) = chebWextrema(Float64, n, k) 
chebWextrema(::Type{T}, n, k) = cospi(T(2*(n-k)+1) / T(2*n+1))

shift_chebWextrema(n, k) = shift_chebWextrema(Float64, n, k) 
shift_chebWextrema(::Type{T}, n, k) where {T} = (chebWextrema(T, n, k) + 1)/2

chebWextrema(n, k) = chebWextrema(Float64, n, k)
function chebWextrema(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    for k = 1:n
        result[k] = shift_chebWextrema(T,n,k)
    end
    return result
end

chebWextrema01(n, k) = chebWextrema01(Float64, n, k)
function chebWextrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    for k = 1:n
        result[k+1] = shift_chebWextrema(T,n,k)
    end
    return result
end

chebWWextrema01(n, k) = chebWWextrema01(Float64, n, k)
function chebWWextrema01(::Type{T}, n) where {T}
    result = Vector{T}(undef, 2*n+2)
    result[1]   = zero(T)
    result[end] = one(T)
    result[2:n+1] = chebWextrema(T,  n, k)
    result[n+2:(2*n+1)] = one(T) .- result[2:n+1]
    sort!(result)
    return result
end
