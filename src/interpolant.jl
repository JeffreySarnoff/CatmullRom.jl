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
within01(n::Int) = collect(0.0:inv(n+1):1.0)[2:end-1]

"""
    uniform01(n::Int)

generate n >= 2 values uniformly spaced from 0..1, includes 0, 1

# Examples
    uniform01(2) == [0.0, 1.0]
    uniform01(3) == [0.0, 0.5, 1.0]
    uniform01(4) == [0.0, 1/3, 2/3, 1.0]
"""
uniform01(n::Int) = collect(LinRange(0.0, 1.0, max(2,n+2)))

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


# julia> pre=CatmullRom.prepoint(xys[1:3]...,);

function prepoint(fn::F, point1::P, point2::P, point3::P) where {F<:Function, P<:OnePoint}
    dimen = length(point1)
    xpre = point1[1] - (point2[1] - point1[1])/8

    point = Vector{typeof(xpre)}(undef, dimen)
    point[1] = xpre
    p1 = point1[1:2]
    p2 = point2[1:2]
    p3 = point3[1:2]
    point[2] = fn(p1, p2, p3, xpre)

    for idx=3:dimen
        p1[2] = point1[idx]
        p2[2] = point2[idx]
        p3[2] = point3[idx]
        point[idx] = fn(p1, p2, p3, xpre)
    end

    return point
end


# julia> post=postpoint(xys[end-2:end]...,);

function postpoint(fn::F, point1::P, point2::P, point3::P) where {F<:Function, P<:OnePoint}
    dimen = length(point1)
    xpost = point3[1] + (point3[1] - point2[1])/8

    point = Vector{typeof(xpost)}(undef, dimen)
    point[1] = xpost
    p1 = point1[1:2]
    p2 = point2[1:2]
    p3 = point3[1:2]
    point[2] = fn(p1, p2, p3, xpost)

    for idx=3:dimen
        p1[2] = point1[idx]
        p2[2] = point2[idx]
        p3[2] = point3[idx]
        point[idx] = fn(p1, p2, p3, xpost)
    end

    return point
end

function linear(pt1, pt2, x)
    a = (pt1[2] - pt2[2]) * x
    b = (pt1[1]*pt2[2] - pt1[2]*pt2[1])
    den = pt1[1] - pt2[1]
    return (a + b) / den
end

function linear(pt1, pt2, pt3, x) # one of the points is ignored
    if (x < pt1[1] < pt2[1]) || (x > pt1[1] > pt2[1])
        res = linear(pt1, pt2, x)
    elseif (pt2[1] < pt3[1] < x) || (pt2[1] > pt3[1] > x)
        res = linear(pt3, pt2, x)
    elseif (pt3[1] < x < pt1[1]) || (pt3[1] > x > pt1[1]) 
        res = linear(pt1, pt3, x)
    else
        throw(DomainError(string(pt1," ",pt2," ",pt3," ",x)))
    end
    return res
end

function quadratic(pt1, pt2, pt3, x)
    t1 = pt2[1] - pt1[1]
    t2 = pt1[1] - pt3[1]
    t3 = pt3[1] - pt2[1]
    t4 = t2 * pt2[2]
    t5 = pt3[1]*pt3[1]
    t6 = pt2[1]*pt2[1]
    t7 = pt1[1]*pt1[1]
    
    s = -(inv(t1) * inv(t2) * inv(t3))    
    
    a = (pt1[2] * t3 + pt3[2] * t1 + t4) * x 
    b = t5 * (pt1[2] - pt2[2])
    c = t6 * (pt3[2] - pt1[2])
    d = t7 * (pt2[2] - pt3[2])
    q = t6 * (pt1[2] * pt3[1] - pt3[2] * pt1[1])
    r = t4 * pt1[1] * pt3[1]
    n = (-pt1[2] * t5 + pt3[2] * t7) * pt2[1]
    aa = a - b - c - d
    bb = r - n - q
    res = (s * (aa * x + bb))

    if any(!isfinite.(res...,))
        res = linear(pt1, pt2, pt3, x)
    end
        
    return res
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
    res = (-(point1[1] - x) * t1 + point1[2])
    
    if any(!isfinite.(res...,))
        res = quadratic(point1, point2, point3, x)
    end
    
    return res
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

