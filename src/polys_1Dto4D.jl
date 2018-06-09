Base.zero(::Type{Poly{T}}) where {T} = Poly(zero(T))
ncoeffs(x::Poly{T}) where {T} = Base.length(x.a)

struct Cubic{T}
    poly::Poly{T}
    δpoly::Poly{T}
        
    function Cubic(c0::T, c1::T, c2::T, c3::T; δx::Bool=false) where {T}
        poly = Poly([c0, c1, c2, c3])
        δpoly = δx ? polyder(poly) : zero(Poly{T})
        return new{T}(poly, δpoly)
    end
end

@inline poly(x::Cubic{T}) where {T} = x.poly
@inline δpoly(x::Cubic{T}) where {T} = ncoeffs(x) === 1 ? polyder(poly(x)) : δpoly(x)
@inline δpoly(x::Cubic{T}, δx::Bool) where {T} = δx ? δpoly(x) : zero(Poly{T})

const dpoly = δpoly

@inline polyval(x::Cubic{T}, z::T) where {T} = polyval(poly(x), z)
@inline δpolyval(x::Cubic{T}, z::T) where {T} = polyval(δpoly(x), z)

const dpolyval = δpolyval


@inline zeros1(::Type{Poly{T}}) where {T} = (zero(Poly{T}),)
@inline zeros2(::Type{Poly{T}}) where {T} = (zero(Poly{T}), zero(Poly{T}))
@inline zeros3(::Type{Poly{T}}) where {T} = (zero(Poly{T}), zero(Poly{T}), zero(Poly{T}))
@inline zeros4(::Type{Poly{T}}) where {T} = (zero(Poly{T}), zero(Poly{T}), zero(Poly{T}), zero(Poly{T}))


struct Cubics{T,N}
    polys::NTuple{N,Poly{T}}
    δpolys::NTuple{N,Poly{T}}
    
    function Cubics{T,N}(x::NTuple{0,T}) where {N, T<:Real}
        polys  = (zero(Poly{T}),)
        δpolys = polys
        return new{T,N}(polys, δpolys)
    end
end

@inline polys(cubics::Cubics{T,N}) where {N,T} = cubics.polys
@inline δpolys(cubics::Cubics{T,N}) where {N,T} = cubics.δpolys

const dpolys = δpolys

@inline polysval(x::Cubics{T,N}, z::T) where {N,T} = polyval.(polys(x), z)
@inline δpolysval(x::Cubics{T,N}, z::T) where {N,T} = polyval.(δpolys(x), z)

@inline function δpolysval(x::Cubics{T,N}, z::T) where {N,T}
    polyders = δpolys(x)
    if any(ncoeffs.(polyders) .== 1)
        polyders = polyder.(polys(x))
    end
    return polyval.(polyders, z)
end

const dpolysval = δpolysval


function Cubics(xcubic::Cubic{T}; δx::Bool=false) where {T}
    polys  = (poly(xcubic),)
    if δx
        δpolys = polyder.(polys)
    else
        δpolys = zeros1(Poly{T})
    end
    return Cubics{T,1}(polys, δpolys)
end

function Cubics(xcubic::Cubic{T}, ycubic::Cubic{T}; δx::Bool=false) where {T}
    polys  = (poly(xcubic), poly(ycubic),)
    if δx
        δpolys = (δpoly(xcubic), δpoly(ycubic),)
    else
        δpolys = zeros2(Poly{T})
    end
    return Cubics{T,2}(polys, δpolys)
end

function Cubics(xcubic::Cubic{T}, ycubic::Cubic{T}, zcubic::Cubic{T}; δx::Bool=false) where {T}
    polys  = (poly(xcubic), poly(ycubic), poly(zcubic))
    if δx
        δpolys = (δpoly(xcubic), δpoly(ycubic), δpoly(zcubic),)
    else
        δpolys = zeros3(Poly{T})
    end
    return Cubics{T,3}(polys, δpolys)
end
        
function Cubics(xcubic::Cubic{T}, ycubic::Cubic{T}, zcubic::Cubic{T}, tcubic::Cubic{T}; δx::Bool=false) where {T}
    polys  = (poly(xcubic), poly(ycubic), poly(zcubic), poly(tcubic))
    if δx
        δpolys = (δpoly(xcubic), δpoly(ycubic), δpoly(zcubic), δpoly(zcubic))
    else
        δpolys = zeros4(Poly{T})
    end
    return Cubics{T,4}(polys, δpolys)
end


"""
    Centripetal Catmull-Rom Cubic polynomials in 1D, 2D, 3D, 4D
    
""" CCR1D, CCR2D, CCR3D, CCR4D

CCR1D = ProtoNT( :x )
CCR2D = ProtoNT( :x, :y )
CCR3D = ProtoNT( :x, :y, :z )
CCR4D = ProtoNT( :x, :y, :z, :t )


# retrieve the Centripetal Catmull-Rom cubic polynomials as Tuples

polys(ccr::CCR1D) = (ccr.x.poly,)
polys(ccr::CCR2D) = (ccr.x.poly, ccr.y.poly)
polys(ccr::CCR3D) = (ccr.x.poly, ccr.y.poly, ccr.z.poly)
polys(ccr::CCR4D) = (ccr.x.poly, ccr.y.poly, ccr.z.poly, ccr.t.poly)

δpolys(ccr::CCR1D) = (ccr.x.δpoly,)
δpolys(ccr::CCR2D) = (ccr.x.δpoly, ccr.y.δpoly)
δpolys(ccr::CCR3D) = (ccr.x.δpoly, ccr.y.δpoly, ccr.z.δpoly)
δpolys(ccr::CCR4D) = (ccr.x.δpoly, ccr.y.δpoly, ccr.z.δpoly, ccr.t.δpoly)

# evaluate Centripetal Catmull-Rom cubic polynomials

for CR in (:CCR1D, :CCR2D, :CCR3D, :CCR4D)
    @eval polyval(ccr::$CR, z::T) where {T} = polyval.(polys(ccr), z)
    @eval δpolyval(ccr::$CR, z::T) where {T} = polyval.(δpolys(ccr), z)
end


#=
 * Compute coefficients for a cubic polynomial
 *   p(s) = c0 + c1*s + c2*s^2 + c3*s^3
 * such that
 *   p(0) = x0, p(1) = x1
 *  and
 *   p'(0) = t0, p'(1) = t1.
=#

function coeffcalc(x0::T, x1::T, t0::T, t1::T) where {T<:Number}
    c0 = x0
    c1 = t0
    c2 = -3*x0 + 3*x1 - 2*t0 - t1
    c3 =  2*x0 - 2*x1 +   t0 + t1
    return Cubic(c0,c1,c2,c3)
end

# compute coefficients for a nonuniform Catmull-Rom spline

function NonuniformCatmullRom(x0::T, x1::T, x2::T, x3::T, dt0::T, dt1::T, dt2::T) where {T}

    # compute tangents when parameterized in [t1,t2]
    t1 = (x1 - x0) / dt0 - (x2 - x0) / (dt0 + dt1) + (x2 - x1) / dt1
    t2 = (x2 - x1) / dt1 - (x3 - x1) / (dt1 + dt2) + (x3 - x2) / dt2

    # rescale tangents for parametrization in [0,1]
    t1 *= dt1
    t2 *= dt1

    return coeffcalc(x1, x2, t1, t2)
end

root4(x) = sqrt(sqrt(x))

function CatmullRomCubics(p0::P, p1::P, p2::P, p3::P) where {T, P<:PT1D{T}}
    d0 = T(root4(Δpoint2(p0, p1)))
    d1 = T(root4(Δpoint2(p1, p2)))
    d2 = T(root4(Δpoint2(p2, p3)))

    # safety check for repeated points
    if (d1 < 1.0e-4)  d1 = T(1.0) end
    if (d0 < 1.0e-4)  d0 = d1 end
    if (d2 < 1.0e-4)  d2 = d1 end

    xcpoly = NonuniformCatmullRom(p0.x, p1.x, p2.x, p3.x, d0, d1, d2)

    return Cubics{T}(xcpoly)
end

function CatmullRomCubics(p0::P, p1::P, p2::P, p3::P) where {T, P<:PT2D{T}}
    d0 = T(root4(Δpoint2(p0, p1)))
    d1 = T(root4(Δpoint2(p1, p2)))
    d2 = T(root4(Δpoint2(p2, p3)))

    # safety check for repeated points
    if (d1 < 1.0e-4)  d1 = T(1.0) end
    if (d0 < 1.0e-4)  d0 = d1 end
    if (d2 < 1.0e-4)  d2 = d1 end

    xcpoly = NonuniformCatmullRom(p0.x, p1.x, p2.x, p3.x, d0, d1, d2)
    ycpoly = NonuniformCatmullRom(p0.y, p1.y, p2.y, p3.y, d0, d1, d2)

    return Cubics{T}(xcpoly, ycpoly)
end

function CatmullRomCubics(p0::P, p1::P, p2::P, p3::P) where {T, P<:PT3D{T}}
    d0 = T(root4(Δpoint2(p0, p1)))
    d1 = T(root4(Δpoint2(p1, p2)))
    d2 = T(root4(Δpoint2(p2, p3)))

    # safety check for repeated points
    if (d1 < 1.0e-4)  d1 = T(1.0) end
    if (d0 < 1.0e-4)  d0 = d1 end
    if (d2 < 1.0e-4)  d2 = d1 end

    xcpoly = NonuniformCatmullRom(p0.x, p1.x, p2.x, p3.x, d0, d1, d2)
    ycpoly = NonuniformCatmullRom(p0.y, p1.y, p2.y, p3.y, d0, d1, d2)
    zcpoly = NonuniformCatmullRom(p0.z, p1.z, p2.z, p3.z, d0, d1, d2)

    return Cubics{T}(xcpoly, ycpoly, zcpoly)
end

function CatmullRomCubics(p0::P, p1::P, p2::P, p3::P) where {T, P<:PT4D{T}}
    d0 = T(root4(Δpoint2(p0, p1)))
    d1 = T(root4(Δpoint2(p1, p2)))
    d2 = T(root4(Δpoint2(p2, p3)))

    # safety check for repeated points
    if (d1 < 1.0e-4)  d1 = T(1.0) end
    if (d0 < 1.0e-4)  d0 = d1 end
    if (d2 < 1.0e-4)  d2 = d1 end

    xcpoly = NonuniformCatmullRom(p0.x, p1.x, p2.x, p3.x, d0, d1, d2)
    ycpoly = NonuniformCatmullRom(p0.y, p1.y, p2.y, p3.y, d0, d1, d2)
    zcpoly = NonuniformCatmullRom(p0.z, p1.z, p2.z, p3.z, d0, d1, d2)
    tcpoly = NonuniformCatmullRom(p0.t, p1.t, p2.t, p3.t, d0, d1, d2)

    return Cubics{T}(xcpoly, ycpoly, zcpoly, tcpoly)
end
