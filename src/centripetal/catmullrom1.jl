struct Cubic{T}
    coeffs::NTuple{4,T}
end
Cubic(a::T, b::T, c::T, d::T) where {T} = Cubic{T}((a,b,c,d))

#=
     Compute coefficients for a cubic polynomial
        p(s) = c0 + c1*s + c2*s^2 + c3*s^3
       such that
          p(0) = x0, p(1) = x1
        and
          p'(0) = t0, p'(1) = t1
=#
function HermiteCubic(x0::T, x1::T, t0::T, t1::T) where {T}
    c0 = x0
    c1 = t0
    c2 = -3 * x0 + 3 * x1 - 2 * t0 - t1
    c3 = 2 * x0 - 2 * x1 + t0 + t1
    return Cubic{T}((c0, c1, c2, c3))
end

# create alt version passing in Vector for results
function centripetal_catmullrom(p0::T, p1::T, p2::T, p3::T; abserr=1.0e-4) where {T}
    dt0 = root4(distancesquared(p0, p1))
    dt1 = root4(distancesquared(p1, p2))
    dt2 = root4(distancesquared(p2, p3))
	
    # safety check for repeated points
    small = eltype(T)(abserr)
    if dt1 < small
        dt1 = one(eltype(T))
    end
    if dt0 < small
        dt0 = dt1
    end
    if dt2 < small
        d2 = dt1
    end
	
    n = length(p0)
    result = Vector{Cubic}(undef, n)
    for i = 1:n
        result[i] = nonuniform_catmullrom(p0[i], p1[i], p2[i], p3[i], dt0, dt1, dt2)
    end
	
    return (result...,)	
end

centripetal_catmullrom(p0::T, p1::T, p2::T, p3::T; abserr=1.0e-4) where {N,T<:NTuple{N,Integer}} =
    centripetal_catmullrom(Float64.(p0), Float64.(p1), Float64.(p2), Float64.(p3); abserr=abserr)

#=
    nonuniform_catmullrom

    compute tangents when parameterized in [t1,t2]
    t1 = ( x1 - x0 ) / dt0 - ( x2 - x0 ) / ( dt0 + dt1 ) + ( x2 - x1 ) / dt1;
    t2 = ( x2 - x1 ) / dt1 - ( x3 - x1 ) / ( dt1 + dt2 ) + ( x3 - x2 ) / dt2;
    Maple simplify
    t1 = ((x2 - x1)*dt0^2 - dt1^2*(-x1 + x0)) / (dt0*(dt0 + dt1)*dt1)
    t2 = ((x3 - x2)*dt1^2 - dt2^2*(-x2 + x1)) / (dt1*(dt1 + dt2)*dt2)
    rescale tangents for parametrization in [0,1]
    t1 *= dt1
    t2 *= dt1
    Maple resimplify
    t1 = ((x2 - x1)*dt0^2 - dt1^2*(-x1 + x0)) / (dt0*(dt0 + dt1))
    t2 = ((x3 - x2)*dt1^2 - dt2^2*(-x2 + x1)) / ((dt1 + dt2)*dt2)
=#
    
function nonuniform_catmullrom(x0::T, x1::T, x2::T, x3::T, dt0::T, dt1::T, dt2::T) where {T}
    dt0sq = dt0 * dt0
    dt1sq = dt1 * dt1
    dt2sq = dt2 * dt2
    
    x2x1 = x2 - x1
    t1den = (dt0 + dt1) * dt0
    t2den = (dt1 + dt2) * dt2
    t1num = (x2x1 * dt0sq) - dt1sq * (x0 - x1)
    t2num = (x2x1 * dt2sq) - dt1sq * (x2 - x3)
    t1 = t1num / t1den
    t2 = t2num / t2den
    
    return HermiteCubic(x1, x2, t1, t2)
end

#=

"""
    catmullrom(points, pointsperarc; extend=true)
    catmullrom(xs, ys, pointsperarc; extend=true)
    catmullrom(xs, ys, zs, pointsperarc; extend=true)

Given abcissa-sequenced path of points (as points or as xs and ys), and
the number of subdivisions to be fit inbetween pair of points
(all neighboring points except the first & last),
obtain the centripetal Catmull-Rom splines that cover each segment.    

By default, extrapolated first and final points
are affixed to the input point sequence so that
all of the input points are present in the resulting sequence
(CatmullRom fitting does not include the first and final points).
If you prefer this not occur, call the function with `extend=false`.

If you prefer to specify the scale factor used in that extrapolation,
use `extendbounds(points, scale=scalefactor)`, and then pass the result
to this function with `extend=false`.
"""
function catmullrom(points::P, pointsperarc::Integer=DefaultPointsPerArc; extend::Bool=true) where P
    catmullrom_requirement(npoints(points))
    if eltype(points) <: Tuple
        apoints = map(x->[x...,], points)
        return catmullrom(apoints, pointsperarc, extend=extend)
    end
    pointsperarc += isodd(pointsperarc)     # force even                                
    
    # ensure that the 'x' values are not coinciding
    changes = norm.(diff(points))[1:end]
    relchanges = changes ./ sum(changes)
    cumrelchanges = cumsum(relchanges)
    pushfirst!(cumrelchanges, 0.0)
    m = vcat(cumrelchanges',reduce(hcat,points))
    m = permutedims(m)
    crpoints = [m[i,:] for i=1:size(m)[1]]

    crpoints = catmullrom_prep(crpoints, pointsperarc, extend=extend)

    # ensure the spline passes through the original values
    return catmullrom_splines(crpoints, pointsperarc)[2:end]
end

function catmullrom_prep(points::P, pointsperarc::Integer=DefaultPointsPerArc; extend::Bool=true) where P
    if isa(points[1], Tuple)
        cr_points = [Float64.([x...,]) for x in points]
    else
        cr_points = [Float64.(x) for x in points]
    end
    if extend
        cr_points = extend_seq(cr_points)
    end
    return cr_points
end

function catmullrom_splines(points::P, pointsperarc::Integer) where P
    coord_type = coordtype(points)
      
    vals_along_each_coord = catmullrom_core(points, pointsperarc) 
    
    if coord_type === Float64
        return vals_along_each_coord
    else
        return [map(coord_type, vals) for vals in vals_along_each_coord]
    end
end

function catmullrom_splines(points::P, pointsperarc::Vector{Integer}) where P
    n_points = npoints(points)
    n_ppa    = length(pointsperarc)
    n_points != n_ppa+1 && throw(ErrorException("length(points) != length(pointsperarc)+1 ($(n_points) != $(n_ppa+1))"))
    
    vals_along_each_coord = []
    for idx = 1:n_ppa
        pts = points[idx:idx+3]
        ppa = pointsperarc[idx]
        vals_each_coord = catmullrom_core(points, ppa) 
        push!(vals_along_each_coord, vals_each_coord)
    end
    return vals_along_each_coord
    
    coord_type = coordtype(points)
    if coord_type === Float64
        return vals_along_each_coord
    else
        return [map(coord_type, vals) for vals in vals_along_each_coord]
    end
end

function catmullrom_core(points::P, pointsperarc::Integer) where P
    catmullrom_requirement(npoints(points))
    n_coords = ncoords(points)
    
    # include the given points (knots) for poly generation
    n_through_points = pointsperarc + 2  # include both endpoints

    #=
        Separate sequences for each coordinate dimension, ordered first to last interpoint span.
        Each poly within a sequence fits an interpoint span for its respective coordinate axis.
        _polys_ is an `(n_points - 3) x (n_coords) 2D array` 
           columns are of one coordinate axis 
           each row holds polys that interpolate one interpoint span across all coordinates.
           (n_points - 2 - 1): drop two extremal points, spans are counted between fenceposts
    =#
    polys = catmullrom_polys(points)
    
    # interpolate using span-relative abcissae [0.0..1.0]
    abcissae01 = range(0.0, stop=1.0, length=n_through_points)

    # vals is an `(n_points - 3) x (n_dims)` array of groups of successive point_along_coord placements
    vals = polyval.(polys, Ref(abcissae01[1:end-1]))
    coord_vals = [collect(Iterators.flatten( vals[:,i] )) for i=1:n_coords]
    endvals = polyval.(polys[end,:], fill(1.0, n_coords))
    @. push!(coord_vals, endvals)
    
    return coord_vals
end

"""
    catmullrom_polys(points)

Traverse a sequence of points (coordinates in 2,3,..n dims),
obtaining interpolatory centripetal Catmull-Rom polynomials.

Yields separate sequences for each coordinate dimension, 
    each sequence is ordered from first to last interpoint span.

Each poly within a sequence fits an interpoint span with respect to
    the sequence's coordinate axis.

_polys_ is an `(n_points - 3) x (n_coords) 2D array`
   columns are of one coordinate axis 
   each row holds polys that interpolate one interpoint span across all coordinates.
"""
function catmullrom_polys(points::P) where P   
    n_points = npoints(points)
    catmullrom_requirement(n_points)

    T = coordtype(points)
                                    # omit the extremal points (-2)
    n_spans  = n_points - 3         # count between fenceposts (-1)
    n_coords = ncoords(points)      # each point has ncoordinates

    polys  = Array{Poly{T}, 2}(undef, n_spans, n_coords)

    for idx = 1:n_spans
        pt₋, pt₀, pt₁, pt₊ = points[idx:idx+3]
        plys = centripetal_catmullrom(pt₋, pt₀, pt₁, pt₊)
        polys[idx,:] = plys
    end

    return polys
end

"""
    centripetal_catmullrom(pt1, pt2, pt3, pt4)

Given four points on a path in their given sequence,
obtain the centripetal Catmull-Rom cubic polynomial
that interpolates the middle two points of the path.
"""
function centripetal_catmullrom(pt₋::T,
                                pt₀::T,
                                pt₁::T,
                                pt₊::T) where T

    tg₀,tg₁     = centripetal_tangents(pt₋, pt₀, pt₁, pt₊)
    c₀,c₁,c₂,c₃ = centripetal_hermite(pt₀, pt₁, tg₀, tg₁)

    coeffs = [[x...] for x in zip(c₀,c₁,c₂,c₃)]
    return Poly.(coeffs)
end

"""
    centripetal_tangents(pt1, pt2, pt3, pt4)

Given four points on a path in their given sequence,
obtain the tangents at the middle two points of the path,
pt2, pt3 delimiting the centripetal curve segment,
each as a unit vector.
"""
function centripetal_tangents(pt₋::T,
                              pt₀::T,
                              pt₁::T,
                              pt₊::T) where T

    # intrasegment centripetal speeds
    dt₋₀ = speed(pt₋, pt₀)
    dt₀₁ = speed(pt₀, pt₁)
    dt₁₊ = speed(pt₁, pt₊)

    dt₋₀, dt₀₁, dt₁₊ = prevent_overlap(T, dt₋₀, dt₀₁, dt₁₊)
    
    # Compute the tangents at the two boundary points
    #   that delimit the interpolatory curve segment.
    #   - tg₀ is the tangent at pt₀
    #   - tg₁ is the tangent at pt₁

    tg₀ = @. (pt₀ - pt₋) /  dt₋₀ -
             (pt₁ - pt₋) / (dt₋₀ + dt₀₁) +
             (pt₁ - pt₀) /  dt₀₁

    tg₁ = @. (pt₁ - pt₀) /  dt₀₁ -
             (pt₊ - pt₀) / (dt₀₁ + dt₁₊) +
             (pt₊ - pt₁) /  dt₁₊

    # Normalize the tangent vectors. These unit vector tangents
    #   provide (guide) centripetal interpolation in [0..1].
    tg₀ = tg₀ .* dt₀₁
    tg₁ = tg₁ .* dt₀₁

    return (tg₀, tg₁)
end

@inline function prevent_overlap(::Type{P}, dt₋₀::T, dt₀₁::T, dt₁₊::T) where {P,T}
    # check if any coordinates coincidee
    ε = coordeps(P)
    # correct for repeated coordinates    
    dt₀₁ = @. ifelse(dt₀₁ <= ε, one(coordtype(P)), dt₀₁)
    dt₋₀ = @. ifelse(dt₋₀ <= ε, dt₀₁, dt₋₀)
    dt₁₊ = @. ifelse(dt₁₊ <= ε, dt₀₁, dt₁₊)
    return dt₋₀, dt₀₁, dt₁₊
 end


"""
    centripetal_hermite(pt₀, pt₁, tg₀, tg₁)
    
- ply(t)  = c0 + c1*t   +   c2*t^2 + c3*t^3    
- dply(t) = c1 + 2*c2*t + 3*c3*t^2    

Given the bounding points `(pt₀, pt₁)` of a  curve segment
and the unit tangent vectors `(tg₀, tg₁)` at those points,
determine the coeffs `(c0, c1, c2, c3)` of cubic `ply(t)`
and of associated quadratic `dply(t)` such that

- ply(0) == pt₀,  dply(0) == tg₀
- ply(1) == pt₁,  dply(1) == tg₁
"""
function centripetal_hermite(pt₀::T, pt₁::T, tg₀::G, tg₁::G) where {T, G}
    c0 = pt₀
    c1 = tg₀
    c2 = @. -3*pt₀ + 3*pt₁ - 2*tg₀ - tg₁
    c3 = @.  2*pt₀ - 2*pt₁ +   tg₀ + tg₁

    return (c0, c1, c2, c3)
end

"""
    speed(pt₀, pt₁)

Obtain the tangent vector along the path
at pt₀ moving to pt₁. The sqrt of that
magnitude is the centripetal speed.
"""
function speed(pt₀::T, pt₁::T) where T
    delta = pt₁ .- pt₀
    return sqrt(sqrtdot(delta, delta))
end
    
@inline function coordeps(::Type{T}) where T   
    coord_type = coordtype(T)
    if coord_type <: AbstractFloat
        ε = eps(coord_type)
    elseif coord_type <: Integer
        ε = zero(coord_type)
    else
        ε = eps(float(dt₀₁))
    end
    return ε
end

sqrtdot(a::T,b::T) where {T} = sqrt(sum(a .* b))

function catmullrom_requirement(n_points::Int)
    n_points >= 4 || throw(DomainError("At least four points are required ($n_points points found)"))
end

=#
