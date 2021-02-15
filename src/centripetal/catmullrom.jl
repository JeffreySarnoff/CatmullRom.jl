const DefaultPointsPerArc = 64

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
    vals = polys.(Ref(abcissae01[1:end-1]))
    coord_vals = [collect(Iterators.flatten( vals[:,i] )) for i=1:n_coords]
    endvals = (polys[end,:]).(fill(1.0, n_coords))
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

    polys  = Array{Polynomial{T}, 2}(undef, n_spans, n_coords)

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
    return Polynomial.(coeffs)
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
