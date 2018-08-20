#=
    fitting functions for simple extrapolation
    one point before the first and
    one point beyond the last
    chosen to provide little curvature mucking
       when used with centripetal CatmullRom
=#

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



const endpointfn = Dict([:Linear => linear,
                         :Quadratic => quadratic,
                         :Thiele3 => thiele3,
                         :Thiele4 => thiele4])
