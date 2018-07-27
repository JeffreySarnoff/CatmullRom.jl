# uniform separation in 0..1 inclusive

"""
    uniformspacing(n::Int)

generate n values uniformly spaced between 0..1 (n does not count including of 0,1) 
"""
uniformspacing(n::Int) = n >= 0 ? collect(0.0:inv((n+2)-1):1.0) : throw(DomainError("$n < 0"))

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


# julia> pre=CatmullRom.prepoint(xys[1:3]...,);

function prepoint(fn::Function, point1, point2, point3)
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

    return (point...,)
end


# julia> post=postpoint(xys[end-2:end]...,);

function postpoint(fn::Function, point1, point2, point3)
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

    return (point...,)
end

function linear(pt1, pt2, pt3, x) # one of the points is ignored
    if x <= pt1[1] 
        res = [pt1 .- (pt1[1] .- x)./(pt2 .- pt1)...,]; res[1] = x
    elseif x >= pt3[1]
        res = [pt3 .+ (x .- pt3[1])./(pt3 .- pt2)...,]; res[1] = x
    else
        throw(DomainError(string(pt1," ",pt2," ",pt3," ",x)))
    end
    return (res...,)
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
        res = linear(point1, point2, point3, x)
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

