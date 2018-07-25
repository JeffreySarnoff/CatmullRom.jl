"""
    catmullrom(points, interpolants)
    `points` is a tuple of points-as-tuples
    `interpolants` is a tuple of values from 0.0 to 1.0 (inclusive)
interpolating points from points[2] through points[end-1] (inclusive)
"""

function catmullrom(points::P, interpolants::Union{NTuple{N,F},Vector{R}}) where
     {I,D,N,R<:Real,F<:Real, P<:Union{NTuple{I, NTuple{D,R}}, Vector{NTuple{D,R}}, Vector{Vector{R}}}}
    npoints = length(points)
    npoints < 4 && throw(ErrorException("at least four points are required"))

    return points === 4 ? catmullrom_4points(points, interpolants) : catmullrom_npoints(points, interpolants)
end

#=
catmullrom(points::Vector{R}, interpolants::NTuple{D,F}) where {D, R, F} =
    catmullrom((points...,), interpolants)

catmullrom(points::NTuple{I, NTuple{N,T}}, interpolants::AbstractArray{F, 1}) where {I, F, N, T} =
    catmullrom(points, (interpolants...,))

catmullrom(points::AbstractArray{T, N}, interpolants::AbstractArray{F, 1}) where {F, N, T} =
    catmullrom((points...,), (interpolants...,))
=#


# ref https://ideone.com/NoEbVM

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


#=
   given four x-ordinate sequenced ND points
   obtain N polys, parameterized over [0,1]
   interpolating from p1 to p2 inclusive
   one poly for each coordinate axis
=#
function catmullrom_polys(points::P) where
     {N,I,D,R<:Real, P<:Union{NTuple{4, NTuple{D,R}}, Vector{NTuple{D,R}}, Vector{Vector{R}}}}
    dt0, dt1, dt2 = prep_centripetal_catmullrom(points)
    pt0, pt1, pt2, pt3 = points

    dimen = length(pt0)
    polys = Vector{Poly{eltype(points[1])}}(undef, dimen)

    for i=1:dimen
        polys[i] = catmullrom_cubic(pt0[i], pt1[i], pt2[i], pt3[i], dt0, dt1, dt2)
    end

    return polys
end

function catmullrom_polys_dpolys(pts::NTuple{4, NTuple{N,T}}) where {N, T}
    polys = catmullrom_polys(pts)
    differentiatepolys = polyder.(polys)
    return polys, differentiatepolys
end

function catmullrom_polys_ipolys(pts::NTuple{4, NTuple{N,T}}) where {N, T}
    polys = catmullrom_polys(pts)
    integratepolys = polyint.(polys)
    return polys, integratepolys
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



qrtrroot(x) = sqrt(sqrt(x))

#=
   determine the delta_traversal constants for the centripetal parameterization
      of the Catmull Rom cubic specified by four points (of increasing abcissae)
=#
function prep_centripetal_catmullrom(points::P) where
     {I,D,R<:Real, P<:Union{NTuple{I, NTuple{D,R}}, Vector{NTuple{D,R}}, Vector{Vector{R}}}}
    dt0 = qrtrroot(dot(points[1], points[2]))
    dt1 = qrtrroot(dot(points[2], points[3]))
    dt2 = qrtrroot(dot(points[3], points[4]))

    #safety check for repeated points
    if (dt1 < 1e-6) dt1 = 1.0 end
    if (dt0 < 1e-6) dt0 = dt1 end
    if (dt2 < 1e-6) dt2 = dt1 end

    return dt0, dt1, dt2
 end


"""
    catmullrom_npoints(points, interpolants)
    `points` is a tuple of points-as-tuples
    `interpolants` is a tuple of values from 0.0 to 1.0 (inclusive)
interpolating points from points[2] through points[end-1] (inclusive)
"""
function catmullrom_npoints(pts::P, interpolants::Union{NTuple{N,F},Vector{F}}) where
     {I,D,N,R,F<:Real, P<:Union{NTuple{I, NTuple{D,R}}, Vector{NTuple{D,R}}, Vector{Vector{R}}}}
    npoints = length(pts)
    points_per_interpolation = length(interpolants)
    totalinterps = (npoints-4+1)*(points_per_interpolation - 1) + 1 # -1 for the shared end|1 point

    dimen = length(pts[1])
    points = Array{R, 2}(undef, (totalinterps,dimen))

    niters = npoints - 5

    points[1:points_per_interpolation,   :] = catmullrom_4points(pts[1:4], interpolants)

    idx₁ = 2; idx₂ = idx₁ + 3; sub₁ = 0; mul₁ = 1
    for k in 1:niters
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
function catmullrom_4points(pts::P, interpolants::Union{NTuple{N,F},Vector{R}}) where
     {D,N,R<:Real,F<:Real, P<:Union{NTuple{4, NTuple{D,R}}, Vector{NTuple{D,R}}, Vector{Vector{R}}}}
    polys = catmullrom_polys(pts)
    ninterps = length(interpolants)
    dimen = length(pts[1])
    points = Array{R, 2}(undef, (ninterps, dimen))
    for col in 1:dimen
        points[1, col] = pts[2][col]
        points[end, col] = pts[3][col]
    end

    ninterps -= 1

    for col in 1:dimen
        ply = polys[col]
        for row in 2:ninterps
            value = interpolants[row]
            points[row, col] = polyval(ply, value)
        end
    end

    return points
end



#=
    Given 4 ND points, roughly approximate the arclength
    of the centripetal Catmull-Rom curvilinear segment
    that would be determined by two bounding points
    and the tangents they determine.

    this algorithm was developed by Jens Gravesen
=#
function approximate_arclength(points::P) where
     {I,D,R<:Real, P<:Union{NTuple{I, NTuple{D,R}}, Vector{NTuple{D,R}}, Vector{Vector{R}}}}
     ldist12 = separation(points[2], points[1])
     ldist23 = separation(points[3], points[2])
     ldist34 = separation(points[4], points[3])
     ldist14 = separation(points[4], points[1])

     linesegments = ldist12 + ldist23 + ldist34
     arclength = (linesegments + ldist14) / 2
     # errorest  = linesegments - ldist14

     return arclength
end

lawofcosines(side1, angle2sides, side2) = side1*side1 + side2*side2 - side1*side2 * 2*cos(angle2sides)

separation(pointa::P, pointb::P) where {N,T, P<:NTuple{N,T}} =
sqrt(lawofcosines(norm(pointa), angle(pointa, pointb), norm(pointb)))

function catmullrom_allpolys(points::P, deriv1::Bool=false, deriv2::Bool=false, integ1::Bool=false) where
     {I,D,R<:Real, P<:Union{NTuple{I, NTuple{D,R}}, Vector{NTuple{D,R}}, Vector{Vector{R}}}}
    polys = catmullrom_polys(points)
    
    d1polys = deriv1 ? polyder.(polys)    : nothing 
    d2polys = deriv2 ? polyder.(d1polys)  : nothing
    i1polys = integ1 ? polyint.(polys)    : nothing

    result = (polys, d1polys, d2polys, i1polys)
    return result
end
