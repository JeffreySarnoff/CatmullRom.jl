__precompile__()

module CentripetalCatmullRom

export catmullrom, into01

using Polynomials
import Polynomials: polyval, polyder

using LinearAlgebra: dot

include("interpolant.jl")


dpolyval(ply::Poly, value::T) where {T} = polyval(polyder(ply), value)


qrtrroot(x) = sqrt(sqrt(x))

# ref https://ideone.com/NoEbVM

#=
   Compute coefficients for a cubic polynomial
     p(s) = c0 + c1*s + c2*s^2 + c3*s^3
   such that
     p(0) = x0, p(1) = x1
    and
     p'(0) = dx0, p'(1) = dx1.
=#
function hermite_cubic(x0::T, x1::T, dx0::T, dx1::T) where {T}
    c0 = x0
    c1 = dx0
    c2 = -3*x0 + 3*x1 - 2*dx0 - dx1
    c3 =  2*x0 - 2*x1 +   dx0 + dx1
    return Poly([c0, c1, c2, c3])
end

# compute nonuniform Catmull-Rom spline over [0, 1]
function catmullrom_cubic(x0::T, x1::T, x2::T, x3::T, dt0::T, dt1::T, dt2::T) where {T}
    # compute tangents when parameterized in [t1,t2]
    t1 = (x1 - x0) / dt0 - (x2 - x0) / (dt0 + dt1) + (x2 - x1) / dt1
    t2 = (x2 - x1) / dt1 - (x3 - x1) / (dt1 + dt2) + (x3 - x2) / dt2
 
    # rescale tangents for parametrization in [0,1]
    t1 *= dt1
    t2 *= dt1
 
    # return hermite cubic over [0,1]
    return hermite_cubic(x1, x2, t1, t2)
end

function prep_centripetal_catmullrom(pts::NTuple{4, NTuple{D,T}}) where {D, T}
    dt0 = qrtrroot(dot(pts[1], pts[2]))
    dt1 = qrtrroot(dot(pts[2], pts[3]))
    dt2 = qrtrroot(dot(pts[3], pts[4]))
 
    #safety check for repeated points
    if (dt1 < 1e-4) dt1 = 1.0 end
    if (dt0 < 1e-4) dt0 = dt1 end
    if (dt2 < 1e-4) dt2 = dt1 end
 
    return dt0, dt1, dt2   
 end

#=
   given four x-ordinate sequenced ND points
   obtain N polys, parameterized over [0,1]
   interpolating from p1 to p2 inclusive
   one poly for each coordinate axis
=#
function catmullrom_polys(pts::NTuple{4, NTuple{N,T}}) where {N, T}
    dt0, dt1, dt2 = prep_centripetal_catmullrom(pts)
    pt0, pt1, pt2, pt3 = pts
    
    polys = Vector{Poly}(undef, N)

    for i=1:N
        polys[i] = catmullrom_cubic(pt0[i], pt1[i], pt2[i], pt3[i], dt0, dt1, dt2)
    end

    return polys
end


function catmullrom_points(pts::NTuple{4, NTuple{D,T}}, interpolants::Union{A,NTuple{N,F}}) where {A<:AbstractArray, N, D, T, F}
    polys = catmullrom_polys(pts)
    ninters = length(interpolants)
    points = Array{T, 2}(undef, (ninters,D))
    for col in 1:D
        ply = polys[col]
        for row in 1:ninters
            value = interpolants[row]
            0.0 <= value <= 1.0 || throw(DomainError("interpolant value ($value) should be in 0.0:1.0"))
            points[row, col] = polyval(ply, value)
        end
    end

    return points
end

catmullrom_points1(pts::NTuple{4, NTuple{D,T}}, interpolants::Union{A,NTuple{N,F}}) where {A<:AbstractArray, N, D, T, F} =
    catmullrom_points(pts, interpolants)[1:end-1,:]

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

#=
    interpolants
    1 2 3 4 5 6         1 2 3 4 5 6
    0         1         1         1
              1 2 3 4 5 6         1 2 3 4 5 6
    1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 6

    (ninterpolants-1)*(npointspans-1) + ninterpolants

    points
    1 2 3 4
      ==    5 
        ==
          ==
 
    p points --> (p-2) interpoland spans, first point and last point are not spanned

    with exactly 4 points, one span of ninterpolants (pt2 and pt3 are included)
    with more points, e.g. 5,  two spans of n-1 then n interpolants

    with k points idx 1:k-3
=#

"""
    catmullrom(points, interpolants)

    `points` is a tuple of points-as-tuples
    `interpolants` is a tuple of values from 0.0 to 1.0 (inclusive)

interpolating points from points[2] through points[end-1] (inclusive)
"""
function catmullrom(points::U1, interpolants::U2) where {U1, U2}
    npoints = length(points)
    npoints < 4 && throw(ErrorException("at least four points are required"))
    
    interps = fixup(interpolants)
    
    if npoints == 4
        return catmullrom_points(points, interps)
    end

    ninterp = length(interpolants)
    ninterp < 2 && throw(ErrorException("at least two interpolants [0,1] are required"))
    ninterp1 = ninterp - 1
    ndimens = length(points[1])
    eltyp = eltype(points[1])

    spans = npoints - 4
    # totalpoints = spans * (ninterp - 1) + ninterp
    totalpoints = (length(points)-1) * (lenth(interpolants)-2) + length(points)


    result = Array{eltyp, ndimens}(undef, (totalpoints, ndimens) )
    
    for i in 0:spans
        k = i+1
        result[(i*ninterp1+1):(k*ninterp1),:] = catmullrom_points1((points[k:k+3]...,), interpolants)
    end
    result[end, :] = [points[end-1]...,]

    return result
end



end # module CentripetalCatmullRom
