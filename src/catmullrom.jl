"""
    catmullrom(points, interpolants)

    `points` is a tuple of points-as-tuples
    `interpolants` is a tuple of values from 0.0 to 1.0 (inclusive)

interpolating points from points[2] through points[end-1] (inclusive)
"""
function catmullrom(points::NTuple{I, NTuple{D,T}}, interpolants::Union{A,NTuple{N,F}}) where {A<:AbstractArray, N, F, I, D, T}
    npoints = length(points)
    npoints < 4 && throw(ErrorException("at least four points are required"))

    return points === 4 ? catmullrom_4points(points, interpolants) : catmullrom_npoints(points, interpolants)
end
 
catmullrom(points::AbstractArray{F, N}, interpolants::NTuple{D,T}) where {F, N, D, T} =
    catmullrom((points...,), interpolants)

catmullrom(points::NTuple{I, NTuple{N,T}}, interpolants::AbstractArray{F, 1}) where {I, F, N, T} =
    catmullrom(points, (interpolants...,))

catmullrom(points::AbstractArray{T, N}, interpolants::AbstractArray{F, 1}) where {F, N, T} =
    catmullrom((points...,), (interpolants...,))



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
function ccr_polys(pts::NTuple{4, NTuple{N,T}}) where {N, T}
    dt0, dt1, dt2 = prep_centripetal_catmullrom(pts)
    pt0, pt1, pt2, pt3 = pts
    
    polys = Vector{Poly{eltype(pts[1])}}(undef, N)

    for i=1:N
        polys[i] = catmullrom_cubic(pt0[i], pt1[i], pt2[i], pt3[i], dt0, dt1, dt2)
    end

    return polys
end

function ccr_polys_dpolys(pts::NTuple{4, NTuple{N,T}}) where {N, T}
    polys = ccr_polys(pts)
    differentiatepolys = polyder.(polys)
    return polys, differentiatepolys
end

function ccr_polys_ipolys(pts::NTuple{4, NTuple{N,T}}) where {N, T}
    polys = ccr_polys(pts)
    integratepolys = polyint.(polys)
    return polys, integratepolys
end

function ccr_allpolys(pts::NTuple{4, NTuple{N,T}}) where {N, T}
    polys = ccr_polys(pts)
    integratepolys = polyint.(polys)
    differentiatepolys = polyder.(polys)"""
    catmullrom(points, interpolants)

    `points` is a tuple of points-as-tuples
    `interpolants` is a tuple of values from 0.0 to 1.0 (inclusive)

interpolating points from points[2] through points[end-1] (inclusive)
"""
function catmullrom(points::NTuple{I, NTuple{D,T}}, interpolants::Union{A,NTuple{N,F}}) where {A<:AbstractArray, N, F, I, D, T}
    npoints = length(points)
    npoints < 4 && throw(ErrorException("at least four points are required"))

    return points === 4 ? catmullrom_4points(points, interpolants) : catmullrom_npoints(points, interpolants)
end
 
catmullrom(points::AbstractArray{F, N}, interpolants::NTuple{D,T}) where {F, N, D, T} =
    catmullrom((points...,), interpolants)

catmullrom(points::NTuple{I, NTuple{N,T}}, interpolants::AbstractArray{F, 1}) where {I, F, N, T} =
    catmullrom(points, (interpolants...,))

catmullrom(points::AbstractArray{T, N}, interpolants::AbstractArray{F, 1}) where {F, N, T} =
    catmullrom((points...,), (interpolants...,))



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
function ccr_polys(pts::NTuple{4, NTuple{N,T}}) where {N, T}
    dt0, dt1, dt2 = prep_centripetal_catmullrom(pts)
    pt0, pt1, pt2, pt3 = pts
    
    polys = Vector{Poly{eltype(pts[1])}}(undef, N)

    for i=1:N
        polys[i] = catmullrom_cubic(pt0[i], pt1[i], pt2[i], pt3[i], dt0, dt1, dt2)
    end

    return polys
end

function ccr_polys_dpolys(pts::NTuple{4, NTuple{N,T}}) where {N, T}
    polys = ccr_polys(pts)
    differentiatepolys = polyder.(polys)
    return polys, differentiatepolys
end

function ccr_polys_ipolys(pts::NTuple{4, NTuple{N,T}}) where {N, T}
    polys = ccr_polys(pts)
    integratepolys = polyint.(polys)
    return polys, integratepolys
end

function ccr_allpolys(pts::NTuple{4, NTuple{N,T}}) where {N, T}
    polys = ccr_polys(pts)
    integratepolys = polyint.(polys)
    differentiatepolys = polyder.(polys)
    return polys, differentiatepolys, integratepolys
end


qrtrroot(x) = sqrt(sqrt(x))

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

 
"""
    catmullrom_npoints(points, interpolants)

    `points` is a tuple of points-as-tuples
    `interpolants` is a tuple of values from 0.0 to 1.0 (inclusive)

interpolating points from points[2] through points[end-1] (inclusive)
"""
function catmullrom_npoints(pts::NTuple{I, NTuple{D,T}}, interpolants::Union{A,NTuple{N,F}}) where {A<:AbstractArray, I, N, D, T, F}
    
    points_per_interpolation = length(interpolants)
    totalinterps = (I-4+1)*(points_per_interpolation - 1) + 1 # -1 for the shared end|1 point
    
    points = Array{T, 2}(undef, (totalinterps,D))
   
    niters = I - 5
   
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
function catmullrom_4points(pts::NTuple{4, NTuple{D,T}}, interpolants::Union{A,NTuple{N,F}}) where {A<:AbstractArray, N, D, T, F}
    polys = ccr_polys(pts)
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


qrtrroot(x) = sqrt(sqrt(x))

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

 
"""
    catmullrom_npoints(points, interpolants)

    `points` is a tuple of points-as-tuples
    `interpolants` is a tuple of values from 0.0 to 1.0 (inclusive)

interpolating points from points[2] through points[end-1] (inclusive)
"""
function catmullrom_npoints(pts::NTuple{I, NTuple{D,T}}, interpolants::Union{A,NTuple{N,F}}) where {A<:AbstractArray, I, N, D, T, F}
    
    points_per_interpolation = length(interpolants)
    totalinterps = (I-4+1)*(points_per_interpolation - 1) + 1 # -1 for the shared end|1 point
    
    points = Array{T, 2}(undef, (totalinterps,D))
   
    niters = I - 5
   
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
   given four ND points in increasing arc-travel order
   and k+2 interpolant values in 0..1 where 0.0 and 1.0 are the `+2`
   determine the k+2 interpolant determined ND points
   where the first interpolant point is the second ND point
   and the final interplant point is the third ND point
=#
function catmullrom_4points(pts::NTuple{4, NTuple{D,T}}, interpolants::Union{A,NTuple{N,F}}) where {A<:AbstractArray, N, D, T, F}
    polys = ccr_polys(pts)
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
    Given 4 ND points, roughly approximate the arclength
    of the centripetal Catmull-Rom curvilinear segment 
    that would be determined by two bounding points
    and the tangents they determine.
    
    this algorithm was developed by Jens Gravesen
=#
function approximate_arclength(pts::NTuple{4, NTuple{N,T}}) where {N, T}
     ldist12 = linearsep(pts[2], pts[1])
     ldist23 = linearsep(pts[3], pts[2])
     ldist34 = linearsep(pts[4], pts[3])
     ldist14 = linearsep(pts[4], pts[1])
     
     linesegments = ldist12 + ldist23 + ldist34
     arclength = (linesegments + ldist14) / 2
     # errorest  = linesegments - ldist14
     
     return arclength   
end

#=
    rough approximation to the length of the arc
       connecting points 2 and 3, this is center arc
       that is interpolated between with each use of
       catmullrom_4points 
=#
function centerarclength(pts::NTuple{4, NTuple{N,T}}) where {N, T}
     ldist12 = linearsep(pts[2], pts[1])
     ldist23 = linearsep(pts[3], pts[2])
     ldist34 = linearsep(pts[4], pts[3])
     ldist14 = linearsep(pts[4], pts[1])
     
     linesegments = ldist12 + ldist23 + ldist34
     arclength = (linesegments + ldist14) / 2
     # errorest  = linesegments - ldist14
     arclength *= (ldist23 / linesegments)
        
     return arclength   
end

# fast, approximate angular separation where both points emanate from the same source (e.g. the origin)
anglesep(pointa, pointb) = acos( dot(pointa, pointb) / sqrt(dot(pointa,pointa) * dot(pointb,pointb)) )
    
    
lawofcosines(side1, angle2sides, side2) = side1*side1 + side2*side2 - side1*side2 * 2*cos(angle2sides)
 
function linearsep(pointa::P, pointb::P) where {N,T, P<:NTuple{N,T}}
    result = sqrt(lawofcosines(norm(pointa), anglesep(pointa, pointb), norm(pointb)))
    isnan(result) ? max(norm(pointa), norm(pointb)) : result
end
   
