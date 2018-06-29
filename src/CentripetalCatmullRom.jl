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

# compute Catmull-Rom cubic curve over [0, 1]
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

#=
   determine the delta_traversal constants for the centripetal parameterization
      of the Catmull Rom cubic specified by four points (of increasing abcissae) 
=#
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

julia> points2D = ([(sin(x*pi),cos(x*pi)) for x=0.0:(0.5/3):0.5]...,)
((0.0, 1.0), (0.49999999999999994, 0.8660254037844387), (0.8660254037844386, 0.5000000000000001), (1.0, 6.123233995736766e-17))

julia> polys = catmullrom_polys(points2D)
2-element Array{Poly,1}:
 Poly(0.49999999999999994 + 0.43301270189221924*x - 0.017949192431122307*x^2 - 0.049038105676658006*x^3)
 Poly(0.8660254037844387 - 0.24999999999999994*x - 0.16506350946109644*x^2 + 0.049038105676658006*x^3)  

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


"""
    catmullrom(points, interpolants)

    `points` is a tuple of points-as-tuples
    `interpolants` is a tuple of values from 0.0 to 1.0 (inclusive)

interpolating points from points[2] through points[end-1] (inclusive)
"""
function catmullrom(points::NTuple{I, NTuple{D,T}}, interpolants::Union{A,NTuple{N,F}}) where {A<:AbstractArray, I, N, D, T, F}
    npoints = length(points)
    npoints < 4 && throw(ErrorException("at least four points are required"))
    
    interps = fixup(interpolants)
    
    if npoints == 4
        return catmullrom_4points(points, interps)
    end
    
   return catmullrom_ipoints(points, interpolants)
end
 
 
"""
    catmullrom_ipoints(points, interpolants)

    `points` is a tuple of points-as-tuples
    `interpolants` is a tuple of values from 0.0 to 1.0 (inclusive)

interpolating points from points[2] through points[end-1] (inclusive)
"""
function catmullrom_ipoints(pts::NTuple{I, NTuple{D,T}}, interpolants::Union{A,NTuple{N,F}}) where {A<:AbstractArray, I, N, D, T, F}
    
    points_per_interpolation = length(interpolants)
    totalinterps = (I-4+1)*(points_per_interpolation - 1) + 1 # -1 for the shared end|1 point
    
    points = Array{T, 2}(undef, (totalinterps,D))
   
    niters = I - 5
   
    points[1:points_per_interpolation,   :] = catmullrom_4points(pts[1:4], interpolants)
   
    idx₁ = 2; idx₂ = idx₁ + 3; sub₁ = 0; mul₁ = 1 
    for k in 1:niters
        global idx₂, idx₁, sub₁, mul₁
        sub₂ = sub₁ - 1
        mul₂ = mul₁ + 1
        points[(mul₁ * points_per_interpolation + sub₁):(mul₂ * points_per_interpolation + sub₂), :] =
            catmullrom_4points(pts[idx₁:idx₂], interpolants)
        idx₁ += 1; idx₂ += 1;
        sub₁, mul₁ = sub₂, mul₂
    end
   
    points[(end-points_per_interpolation+1):end, :] = catmullrom_4points(pts[end-3:end], interpolants)
   
    return points
end

#=
   given four ND points in increasing abcissae order
   and k+2 interpolant values in 0..1 where 0.0 and 1.0 are the +2
   determine the k+2 interpolant determined ND points
   where the first interpolant point is the second ND point
   and the final interplant point is the third ND point
=#
function catmullrom_4points(pts::NTuple{4, NTuple{D,T}}, interpolants::Union{A,NTuple{N,F}}) where {A<:AbstractArray, N, D, T, F}
    validate_interpolants(interpolants)
    
    polys = catmullrom_polys(pts)
    ninterps = length(interpolants)
    
    points = Array{T, 2}(undef, (ninterps,D))
    for col in 1:D
        points[1, col] = pts[2][col]
        points[end, col] = pts[3][col]
    end
   
    ninterps -= 1
   
    for col in 1:D
        ply = polys[col]
        for row in 2:ninterps
            value = interpolants[row]
            points[row, col] = polyval(ply, value)
        end
    end

    return points
end

#=

function catmullrom_5points(pts::NTuple{5, NTuple{D,T}}, interpolants::Union{A,NTuple{N,F}}) where {A<:AbstractArray, N, D, T, F}
    
    points_per_interpolation = length(interpolants)
    totalinterps = (5-4+1)*(points_per_interpolation - 1) + 1 # -1 for the shared end|1 point
   
    points = Array{T, 2}(undef, (totalinterps,D))
   
    points[1:points_per_interpolation,   :] = catmullrom_4points(pts[1:4], interpolants)
    points[(end-points_per_interpolation+1):end, :] = catmullrom_4points(pts[2:5], interpolants)
   
    return points
end

function catmullrom_6points(pts::NTuple{6, NTuple{D,T}}, interpolants::Union{A,NTuple{N,F}}) where {A<:AbstractArray, N, D, T, F}
    
    points_per_interpolation = length(interpolants)
    totalinterps = (6-4+1)*(points_per_interpolation - 1) + 1 # -1 for the shared end|1 point
   
    points = Array{T, 2}(undef, (totalinterps,D))
   
    points[1:points_per_interpolation,   :] = catmullrom_4points(pts[1:4], interpolants)
    points[points_per_interpolation:(2*points_per_interpolation-1),   :] = catmullrom_4points(pts[2:5], interpolants)
    points[(end-points_per_interpolation+1):end, :] = catmullrom_4points(pts[3:6], interpolants)
   
    return points
end

function catmullrom_7points(pts::NTuple{7, NTuple{D,T}}, interpolants::Union{A,NTuple{N,F}}) where {A<:AbstractArray, N, D, T, F}
    
    points_per_interpolation = length(interpolants)
    totalinterps = (7-4+1)*(points_per_interpolation - 1) + 1 # -1 for the shared end|1 point
   
    points = Array{T, 2}(undef, (totalinterps,D))
   
    points[1:points_per_interpolation,   :] = catmullrom_4points(pts[1:4], interpolants)
    points[(1*points_per_interpolation):(2*points_per_interpolation-1),   :] = catmullrom_4points(pts[2:5], interpolants)
    points[(2*points_per_interpolation-1):(3*points_per_interpolation-2),   :] = catmullrom_4points(pts[3:6], interpolants)
    points[(end-points_per_interpolation+1):end, :] = catmullrom_4points(pts[4:7], interpolants)
   
    return points
end

function catmullrom_8points(pts::NTuple{8, NTuple{D,T}}, interpolants::Union{A,NTuple{N,F}}) where {A<:AbstractArray, N, D, T, F}
    
    points_per_interpolation = length(interpolants)
    totalinterps = (8-4+1)*(points_per_interpolation - 1) + 1 # -1 for the shared end|1 point
   
    points = Array{T, 2}(undef, (totalinterps,D))
   
    points[1:points_per_interpolation,   :] = catmullrom_4points(pts[1:4], interpolants)
    points[(1*points_per_interpolation-0):(2*points_per_interpolation-1),   :] = catmullrom_4points(pts[2:5], interpolants)
    points[(2*points_per_interpolation-1):(3*points_per_interpolation-2),   :] = catmullrom_4points(pts[3:6], interpolants)
    points[(3*points_per_interpolation-2):(4*points_per_interpolation-3),   :] = catmullrom_4points(pts[4:7], interpolants)
    points[(end-points_per_interpolation+1):end, :] = catmullrom_4points(pts[5:8], interpolants)
   
    return points
end

=#
   
#=
    npoints = length(points)
    npoints < 4 && throw(ErrorException("at least four points are required"))
    
    interps = fixup(interpolants)
    
    if npoints == 4
        return catmullrom_4points(points, interps)
    end
=#

#=
# all but the last from catmullrom_points (all except the third ND point)
catmullrom_4points1(pts::NTuple{4, NTuple{D,T}}, interpolants::Union{A,NTuple{N,F}}) where {A<:AbstractArray, N, D, T, F} =
    catmullrom_4points(pts, interpolants)[1:end-1,:]


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
        return catmullrom_4points(points, interps)
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
=#

end # module CentripetalCatmullRom
