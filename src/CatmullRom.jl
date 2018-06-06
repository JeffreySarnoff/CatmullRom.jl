# q.v. https://stackoverflow.com/questions/9489736/catmull-rom-curve-with-no-cusps-and-no-self-intersections/23980479#23980479

module CatmullRom

export CentripetalCatmullRom, Point2D

using Polynomials
import Polynomials: Poly, polyval, polyint, polyder

struct Point2D{T}
    x::T
    y::T
end

struct Point3D{T}
    x::T
    y::T
    z::T
end

function DistSquared(p::Point2D, q::Point2D)
    dx = q.x - p.x
    dy = q.y - p.y
    return dx*dx + dy*dy
end

function DistSquared(p::Point3D, q::Point3D)
    dx = q.x - p.x
    dy = q.y - p.y
    dz = q.z - p.z
    return dx*dx + dy*dy + dz*dz
end

struct CubicPoly{T}
    c0::T
    c1::T
    c2::T
    c3::T
    
    poly::Poly{T}
    diff1poly::Poly{T}
    diff2poly::Poly{T}
    integratedpoly::Poly{T}
    
    function CubicPoly(c0::T, c1::T, c2::T, c3::T) where {T}
        poly = Poly([c0, c1, c2, c3])
        diff1poly = polyder(poly)
        diff2poly = polyder(diff1poly)
        integratedpoly = polyint(poly)
        return new{T}(c0, c1, c2, c3, poly, diff1poly, diff2poly, integratedpoly)
    end
    
end

@inline poly(x::CubicPoly{T}) where {T} = x.poly
@inline diff1poly(x::CubicPoly{T}) where {T} = x.diff1poly
@inline diff2poly(x::CubicPoly{T}) where {T} = x.diff2poly
@inline integratedpoly(x::CubicPoly{T}) where {T} = x.integratedpoly

@inline polyval(x::CubicPoly{T}, z::T) where {T} = polyval(x.poly, z)
@inline diff1polyval(x::CubicPoly{T}, z::T) where {T} = polyval(x.diff1poly, z)
@inline diff2polyval(x::CubicPoly{T}, z::T) where {T} = polyval(x.diff2poly, z)
@inline integratedpolyval(x::CubicPoly{T}, z::T) where {T} = polyval(x.integratedpoly, z)

function evalpoly(cpoly::CubicPoly{T}, t::T) where {T<:Number}
    t1 = t
    t2 = t1 * t1
    t3 = t2 * t1
    t3 *= cpoly.c3
    t2 *= cpoly.c2
    t1 *= cpoly.c1
    t1 += cpoly.c0
    t1 += t2
    t1 += t3
    return t1
end

evalpoly(cpoly::CubicPoly{Float32}, t::Float64) = polyval(cpoly, Float32(t))

#=
 * Compute coefficients for a cubic polynomial
 *   p(s) = c0 + c1*s + c2*s^2 + c3*s^3
 * such that
 *   p(0) = x0, p(1) = x1
 *  and
 *   p'(0) = t0, p'(1) = t1.
=#

function InitCubicPoly(x0::T, x1::T, t0::T, t1::T) where {T<:Number}
    c0 = x0
    c1 = t0
    c2 = -3*x0 + 3*x1 - 2*t0 - t1
    c3 =  2*x0 - 2*x1 +   t0 + t1
    return CubicPoly(c0,c1,c2,c3)
end

# compute coefficients for a nonuniform Catmull-Rom spline
function NonuniformCatmullRom(x0::T, x1::T, x2::T, x3::T, dt0::T, dt1::T, dt2::T) where {T}
    # compute tangents when parameterized in [t1,t2]
    t1 = (x1 - x0) / dt0 - (x2 - x0) / (dt0 + dt1) + (x2 - x1) / dt1
    t2 = (x2 - x1) / dt1 - (x3 - x1) / (dt1 + dt2) + (x3 - x2) / dt2

    # rescale tangents for parametrization in [0,1]
    t1 *= dt1
    t2 *= dt1

    return InitCubicPoly(x1, x2, t1, t2)
end



function CentripetalCatmullRom(p0::Point2D{T}, p1::Point2D{T}, p2::Point2D{T}, p3::Point2D{T}) where {T}
    dt0 = T(sqrt(sqrt(DistSquared(p0, p1))))
    dt1 = T(sqrt(sqrt(DistSquared(p1, p2))))
    dt2 = T(sqrt(sqrt(DistSquared(p2, p3))))

    # safety check for repeated points
    if (dt1 < 1.0e-4)  dt1 = T(1.0) end
    if (dt0 < 1.0e-4)  dt0 = dt1 end
    if (dt2 < 1.0e-4)  dt2 = dt1 end

    xcpoly = NonuniformCatmullRom(p0.x, p1.x, p2.x, p3.x, dt0, dt1, dt2)
    ycpoly = NonuniformCatmullRom(p0.y, p1.y, p2.y, p3.y, dt0, dt1, dt2)

    return xcpoly, ycpoly
end


function CentripetalCatmullRom(p0::Point3D{T}, p1::Point3D{T}, p2::Point3D{T}, p3::Point3D{T}) where {T}
    dt0 = T(sqrt(sqrt(DistSquared(p0, p1))))
    dt1 = T(sqrt(sqrt(DistSquared(p1, p2))))
    dt2 = T(sqrt(sqrt(DistSquared(p2, p3))))

    # safety check for repeated points
    if (dt1 < 1.0e-4)  dt1 = T(1.0) end
    if (dt0 < 1.0e-4)  dt0 = dt1 end
    if (dt2 < 1.0e-4)  dt2 = dt1 end

    xcpoly = NonuniformCatmullRom(p0.x, p1.x, p2.x, p3.x, dt0, dt1, dt2)
    ycpoly = NonuniformCatmullRom(p0.y, p1.y, p2.y, p3.y, dt0, dt1, dt2)
    zcpoly = NonuniformCatmullRom(p0.z, p1.z, p2.z, p3.z, dt0, dt1, dt2)

    return xcpoly, ycpoly, zpoly
end

function test()
    p0 = Point2D(0.0, 0.0)
    p1 = Point2D(1.0, 1.0)
    p2 = Point2D(1.5, 1.25)
    p3 = Point2D(2.0, 0.0)

    xcpoly, ycpoly = CentripetalCatmullRom(p0, p1, p2, p3)

    for i=0:10
        xcoord = evalpoly(xcpoly, 0.1*i)
        ycoord = evalpoly(ycpoly, 0.1*i)
        coord = Point2D(xcoord, ycoord)
        println(coord)
    end
end


end # module
