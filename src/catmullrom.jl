"""
    catmullrom(points, interpolants; endpoints::Symbol=Thiele3)
    catmullrom(points, ninterpolants::Int; endpoints::Symbol=Thiele3)

    `points` is a tuple/vector of points-as-tuples/vectors
    `interpolants` is a tuple/vector of values from 0.0 to 1.0 (inclusive)

all of the interpolants are applied to each segment

endpoints==Omit: interpolating from points[2] through points[end-1]
endpoints==Thiele3:  interpolating from points[1] through points[end]
"""
function catmullrom(points::PointSeq, interpolants::ValueSeq; allpoints::Bool=true) where {M,D,R,L,F}
    npoints = length(points)
    npoints < 4 && throw(ErrorException("at least four points are required"))

    return if endpoints === Thiele3
               T = typeof(points[1]) <: Tuple ? Tuple : Vector
               catmullrom_npoints(augmentends(T, points), interpolants)
           elseif npoints > 4
               catmullrom_npoints(points, interpolants)
           else
               catmullrom_4points(points, interpolants)
           end
end

function catmullrom(points::PointSeq, ninterpolants::Int; allpoints::Bool=true) where {M,D,R}
    interpolants = uniformspacing(ninterpolants)
    return catmullrom(points, interpolants, allpoints)
end

@inline function augmentends(::Type{Tuple}, points::PointSeq) where {M,D,R}
    pre  = prepoint(points[1:3]...,)
    post = postpoint(points[end-2:end]...,)

    return [pre, points..., post]
end

@inline function augmentends(::Type{Vector}, points::PointSeq) where {M,D,R}
    pre  = [prepoint(points[1:3]...,)...,]
    post = [postpoint(points[end-2:end]...,)...,]

    return [pre, points..., post]
end


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

qrtrroot(x) = sqrt(sqrt(x))

#=
   determine the delta_traversal constants for the centripetal parameterization
      of the Catmull Rom cubic specified by four points (of increasing abcissae)
=#
function prep_centripetal_catmullrom(points::PointSeq) where {M,D,R}
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
function catmullrom_npoints(pts::PointSeq, interpolants::ValueSeq) where {M,D,R,L,F}
    npoints = length(pts)
    points_per_interpolation = length(interpolants)
    totalinterps = (npoints-4+1)*(points_per_interpolation - 1) + 1 # -1 for the shared end|1 point

    dimen = length(pts[1])
    T = NTuple{dimen, typeof(pts[1][1])}
    points = Array{T, 1}(undef, totalinterps)

    niters = npoints - 5

    points[1:points_per_interpolation, :] .= catmullrom_4points(pts[1:4], interpolants)

    idx₁ = 2; idx₂ = idx₁ + 3; sub₁ = 0; mul₁ = 1
    for k in 1:niters
        sub₂ = sub₁ - 1
        mul₂ = mul₁ + 1
        points[(mul₁ * points_per_interpolation + sub₁):(mul₂ * points_per_interpolation + sub₂), :] .=
            catmullrom_4points(pts[idx₁:idx₂], interpolants)
        idx₁ += 1; idx₂ += 1;
        sub₁, mul₁ = sub₂, mul₂
    end

    points[(end-points_per_interpolation+1):end, :] .= catmullrom_4points(pts[end-3:end], interpolants)

    return points
end

#=
   given four ND points in increasing abcissae order
   and k+2 interpolant values in 0..1 where 0.0 and 1.0 are the +2
   determine the k+2 interpolant determined ND points
   where the first interpolant point is the second ND point
   and the final interplant point is the third ND point
=#
function catmullrom_4points(pts::PointSeq, interpolants::ValueSeq) where {M,D,R,L,F}
    polys = catmullrom_polys(pts)
    totalinterps = length(interpolants)
    dimen = length(pts[1])
    T = typeof(pts[1][1])
    points = Array{T, 2}(undef, (totalinterps,dimen))

    for col in 1:dimen
        points[1, col] = pts[2][col]
        points[end, col] = pts[3][col]
    end

    totalinterps -= 1

    for col in 1:dimen
        ply = polys[col]
        for row in 2:totalinterps
            value = interpolants[row]
            points[row, col] = polyval(ply, value)
        end
    end

    return [(points[i,:]...,) for i=1:size(points)[1]]
end

#=
   given four x-ordinate sequenced ND points
   obtain N polys, parameterized over [0,1]
   interpolating from p1 to p2 inclusive
   one poly for each coordinate axis
=#
function catmullrom_polys(points::PointSeq) where {M,D,R}
    length(points) == 4 || throw(DomainError("exactly four points are required"))

    dt0, dt1, dt2 = prep_centripetal_catmullrom(points)
    pt0, pt1, pt2, pt3 = points

    dimen = length(pt0)
    polys = Vector{Poly{eltype(points[1])}}(undef, dimen)

    for i=1:dimen
        polys[i] = catmullrom_cubic(pt0[i], pt1[i], pt2[i], pt3[i], dt0, dt1, dt2)
    end

    return polys
end

function catmullrom_allpolys(points::PointSeq, deriv1::Bool=false, deriv2::Bool=false, integ1::Bool=false) where {M,D,R}
    polys = catmullrom_polys(points)

    d1polys = deriv1 ? polyder.(polys)    : nothing
    d2polys = deriv2 ? polyder.(d1polys)  : nothing
    i1polys = integ1 ? polyint.(polys)    : nothing

    result = (polys, d1polys, d2polys, i1polys)
    return result
end
