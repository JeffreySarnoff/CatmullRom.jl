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

function prepoint(point1, point2, point3)
    dimen = length(point1)
    xpre = point1[1] - (point2[1] - point1[1])/8
    
    point = Vector{typeof(xpre)}(undef, dimen)
    point[1] = xpre
    p1 = point1[1:2]
    p2 = point2[1:2]
    p3 = point3[1:2]
    point[2] = thiele3(p1, p2, p3, xpre)
    
    for idx=3:dimen
        p1[2] = point1[idx]
        p2[2] = point2[idx]
        p3[2] = point3[idx]
        point[idx] = thiele3(p1, p2, p3, xpre)
    end
    
    return point
end


function postpoint(point1, point2, point3)
    dimen = length(point1)
    xpost = point3[1] + (point3[1] - point2[1])/8
    
    point = Vector{typeof(xpre)}(undef, dimen)
    point[1] = xpost
    p1 = point1[1:2]
    p2 = point2[1:2]
    p3 = point3[1:2]
    point[2] = thiele3(p1, p2, p3, xpost)
    
    for idx=3:dimen
        p1[2] = point1[idx]
        p2[2] = point2[idx]
        p3[2] = point3[idx]
        point[idx] = thiele3(p1, p2, p3, xpost)
    end
    
    return point
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
    t2 = inv(t2)
    t1 = -(point2[1] - x) * t2 + t1
    t1 = inv(t1)
    return (-(point1[1] - x) * t1 + point1[2])
end

function thiele4(point1, point2, point3, point4, x)
    t1 = point1[2] - point2[2]
    t2 = point2[2] - point3[2]
    t2 = inv(t2)
    t1 = inv(t1)
    t2 = (point2[1] - point3[1]) * t2
    t1 = (point1[1] - point2[1]) * t1
    t3 = -t2 + t1
    t4 = point3[2] - point4[2]
    t4 = inv(t4)
    t4 = -(point3[1] - point4[1]) * t4 + t2
    t4 = inv(t4)
    t3 = inv(t3)
    t3 = (point1[1] - point3[1]) * t3
    t4 = -(point2[1] - point4[1]) * t4 + point2[2] - point3[2] + t3
    t4 = inv(t4)
    t2 = -(point1[1] - point4[1]) * t4 + t1 - t2
    t2 = inv(t2)
    t2 = (point3[1] - x) * t2 - point1[2] + point2[2] + t3
    t2 = inv(t2)
    t1 = -(point2[1] - x) * t2 + t1
    t1 = inv(t1)
    return (-(point1[1] - x) * t1 + point1[2])
end

# uniform separation in 0..1 inclusive
uniformsep(n::Int) = n >= 2 ? collect(0.0:inv(n-1):1.0) : throw(DomainError("$n < 2"))
