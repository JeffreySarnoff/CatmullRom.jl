#=
    Given a centripetal CatmullRom with
    its 4 spline control points (CR[1:4]),
    where the spline segment is a cubic arc
    from CR[2] to CR[3] where the arc is
    interpolatory at CR[2] and CR[3],
    obtain
    the equivalent spline segment as a
    Bezier cubic arc with its 4 spline
    control points (BZ[1:4]) where the arc is
    noninterpolatory and within the convex hull
    formed of the three line segments
    BZ[1]..BZ[2], BZ[2]..BZ[3], BZ[3]..BZ[4].
=#

abstract type AbstractCubicPoly end
abstract type AbstractCubicBezier <: AbstractCubicPoly end
abstract type AbstractCatmullRom <: AbstractCubicBezier end

struct CubicPoly{T} <: AbstractCubicPoly
    c0::T
    c1::T
    c2::T
    c3::T

    CubicPoly(c::NTuple{4,T}) where {T} =
        new{T}(c...)
end

struct CatmullRom{T} <: AbstractCatmullRom
    pA::T
    p1::T
    p2::T
    pZ::T

    CatmullRom(p::NTuple{4,T}) where {T} =
        new{T}(p...)
end

struct CubicBezier{T} <: AbstractCubicBezier
    b0::T
    b1::T
    b2::T
    b3::T

    CubicBezier(b::NTuple{4,T}) where {T} =
        new{T}(b...)
end

# Given a CatmullRom [cubic] arc,
# obtain a BezierCubic arc
#
#   this version presupposes the centripetal
#   parameterization of a CatmullRom cubic,
#   so alpha is implicit and alpha == 1/2
#
# from ref [1]:
# Let p0,p1,p2,p3 be the four consecutive control points
# of the Catmull-Rom curve with corresponding consecutive
# four parameter values v0, v1, v2, v3:
#   0, d01α, d12α + d01α, d23α + d12α + d01α
#   where dij = norm(pj .- pi), dijα = dij^α
#
# The control points of the cubic Bezier curve Bj
# (j ∈ {0, 1, 2, 3}) representing this polynomial
# between d01^α and d12^α + d01^α,
# reparameterized to lie in the range [0,1]
# are calculated as b0, b1, b2, b3.
#
# ref:
#
# [1] On the Parameterization of Catmull-Rom Curves
# by Cem Yuksel, Scott Schaefer, John Keyser
# 2009 SIAM/ACM Joint Conference on Geometric and Physical Modeling
# http://www.cemyuksel.com/research/catmullrom_param/catmullrom.pdf
#
# [2] Paramterizations and Applications of Catmull-Rom Curves
# by Cem Yuksel, Scott Schaefer, John Keyser
# this code derived from equations on pg6
# http://www.cemyuksel.com/research/catmullrom_param/catmullrom_cad.pdf
#
function BezierCubic(x::CatmullRom{T}) where {T}

    # p0,p1,p2,p3 are the CatmullRom control points
    #
    p0 = x.cA
    p1 = x.c1
    p2 = x.c2
    p3 = x.cZ

    # d01, d12, d23 are distances between
    # adjacent [CatmullRom control] points
    # (p0, p1), (p1, p2), (p2, p3)
    #
    d01 = norm(p1 .- p0)
    d12 = norm(p2 .- p1)
    d23 = norm(p3 .- p2)

    # d01α, d12α, d23α are
    # distance^alpha with alpha==1/2
    # and distances between ajacent
    # CatmullRom control points
    # (p0, p1), (p1, p2), (p2, p3)
    #
    d01α = sqrt(d01)
    d12α = sqrt(d12)
    d23α = sqrt(d23)

    # b0, b1, b2, b3 are cubic Bezier control points
    #
    b0 = p1
    b3 = p2
    # b1den = 3*d01α * (d01α + d12α)
    # b2den = 3*d23α * (d23α + d12α)
    b1den = 3*(d01α*d01α + d01α*d12α)
    b2den = 3*(d23α*d23α + d23α*d12α)
    b1num = (d01 .* p2) - (d12 .* p0) + ((2*d01 + 3*d01α*d12α + d12) .* p1)
    b2num = (d23 .* p1) - (d12 .* p3) + ((2*d23 + 3*d23α*d12α + d12) .* p2)
    b1 = b1num / b1den
    b2 = b2num / b2den

    return BezierCubic((b0, b1, b2, b3))
end

# Given a CatmullRom [cubic] arc,
# obtain a BezierCubic arc
#
#   this version presupposes the centripital
#   parameterization of a CatmullRom cubic,
#   so alpha is implicit and alpha == 1/2
#
# from ref [1]:
# Let p0,p1,p2,p3 be the four consecutive control points
# of the Catmull-Rom curve with corresponding consecutive
# four parameter values v0, v1, v2, v3:
#   0, d01α, d12α + d01α, d23α + d12α + d01α
#   where dij = norm(pj .- pi), dijα = dij^α
#
# The control points of the cubic Bezier curve Bj
# (j ∈ {0, 1, 2, 3}) representing this polynomial
# between d01^α and d12^α + d01^α,
# reparameterized to lie in the range [0,1]
# are calculated as b0, b1, b2, b3.
#
# ref:
#
# [1] On the Parameterization of Catmull-Rom Curves
# by Cem Yuksel, Scott Schaefer, John Keyser
# 2009 SIAM/ACM Joint Conference on Geometric and Physical Modeling
# http://www.cemyuksel.com/research/catmullrom_param/catmullrom.pdf
#
# [2] Paramterizations and Applications of Catmull-Rom Curves
# by Cem Yuksel, Scott Schaefer, John Keyser
# this code derived from equations on pg6
# http://www.cemyuksel.com/research/catmullrom_param/catmullrom_cad.pdf
#
function BezierCubic(x::CatmullRom{T}) where {T}

    # p0,p1,p2,p3 are the CatmullRom control points
    #
    p0 = x.cA
    p1 = x.c1
    p2 = x.c2
    p3 = x.cZ

    # d01, d12, d23 are distances between
    # adjacent [CatmullRom control] points
    # (p0, p1), (p1, p2), (p2, p3)
    #
    d01 = norm(p1 .- p0)
    d12 = norm(p2 .- p1)
    d23 = norm(p3 .- p2)

    # d01α, d12α, d23α are
    # distance^alpha with alpha==1/2
    # and distances between ajacent
    # CatmullRom control points
    # (p0, p1), (p1, p2), (p2, p3)
    #
    d01α = sqrt(d01)
    d12α = sqrt(d12)
    d23α = sqrt(d23)

    # b0, b1, b2, b3 are cubic Bezier control points
    #
    b0 = p1
    b3 = p2
    # b1den = 3*d01α * (d01α + d12α)
    # b2den = 3*d23α * (d23α + d12α)
    b1den = 3*(d01α*d01α + d01α*d12α)
    b2den = 3*(d23α*d23α + d23α*d12α)
    b1num = (d01 .* p2) - (d12 .* p0) + ((2*d01 + 3*d01α*d12α + d12) .* p1)
    b2num = (d23 .* p1) - (d12 .* p3) + ((2*d23 + 3*d23α*d12α + d12) .* p2)
    b1 = b1num / b1den
    b2 = b2num / b2den

    return BezierCubic((b0, b1, b2, b3))
end
