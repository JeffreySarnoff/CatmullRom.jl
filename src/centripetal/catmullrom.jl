"""
    catmullrom(n_more_points, points_given; extend=true)
    catmullrom(n_more_points, xs, ys; extend=true)
    catmullrom(n_more_points, xs, ys, zs; extend=true)

Given abcissa-sequenced path of points (as points or as xs and ys), and
the number of subdivisions to be fit inbetween
each path-adjacent, non-extremal pair of points
(all neighboring points except the first & last),
obtain the centripetal Catmull-Rom splines that
cover each interpoint segment between neighbors.    

Use these cubic polynomials to obtain coordinates
for each bounding point of an interpoint segment,
and provide them with the given points to obtain
an augmented abcissa-sequenced path of points.

By default, extrapolated first and final points
are affixed to the input point sequence so that
all of the input points are present in the resulting sequence
(CatmullRom fitting does not include the first and final points).
If you prefer this not occur, call the function with
`extend=false`.

If you prefer to specify the scale factor used
in that extrapolation, use `extendbounds(points, scale=scalefactor)`,
and then pass that result to this function with `extend=false`.
"""
function catmullrom(n_more_points::Integer, points::Points; extend::Bool=true)
    catmullrom_requirement(npoints(points))    
    n_interpolants = n_interpolants + isodd(n_interpolants)  # force even
    
    if points[1] == points[end]                              # curve is closed
        pushfirst!(push!(points, points[2]), points[end-1])  #   close the spline
    elseif extend                                            # curve is open 
        points = extendbounds(points)                        #   cap the spline 
    end    
    
    return catmullrom_splines(n_more_points, points)
end

function catmullrom(n_more_points::Integer, xs::Vector{T}, ordinates::Vararg{Vector{T}}; extend::Bool=true) where {T<:Points}
    n_points = npoints(xs)
    catmullrom_requirement(n_points)
    all(n_points .== length.(ordinates)) || throw(DomainError("lengths must match ($n_points, $(length.(ordinates)))"))
    
    n_more_points = max(n_more_points, 2*(n_points -1))
    return catmullrom(n_more_points, zip(xs, ordinates...), extend=extend)
end

    
catmullrom(n_more_points::Integer, points::Base.Iterators.Zip; extend::Bool=true) =
    catmullrom(n_more_points, collect(points), extend=extend)


function catmullrom_splines(n_more_points::Integer, points::Points)
    coord_type = coordtype(points)
      
    vals_along_each_coord = catmullrom_core(n_more_points, points) 
    
    if coord_type === Float64
        return vals_along_each_coord
    else
        return [map(coord_type, vals) for vals in vals_along_each_coord]
    end
end


function catmullrom_core(n_more_points::Integer, points::Points)
    catmullrom_requirement(points)
    
    n_points = npoints(points)
    n_coords = ncoords(points)
    
    # include the given points (knots) for poly generation
    n_through_points = n_more_points + 2  # include both endpoints

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
    abcissae01 = range(0.0, 1.0, length=n_through_points)

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
function catmullrom_polys(points::Points)
    catmullrom_requirement(points)
    
    n_points = npoints(points)
    
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
                                pt₊::T) where {T<:OnePoint}

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
                              pt₊::T) where {T<:OnePoint}

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

@inline function prevent_overlap(::Type{P}, dt₋₀::T, dt₀₁::T, dt₁₊::T) where {T, P<:OnePoint}
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
function centripetal_hermite(pt₀::T, pt₁::T, tg₀::G, tg₁::G) where {T<:OnePoint, G}
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
function speed(pt₀::T, pt₁::T) where {T<:OnePoint}
    delta = pt₁ .- pt₀
    return sqrt(sqrtdot(delta, delta))
end
    
@inline function coordeps(::Type{T}) where {T<:OnePoint}    
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
