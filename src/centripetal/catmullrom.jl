"""
    catmullrom(  points_along_curve, n_between_points )
               ; augment = true, iterator = false   )

Given abcissa-sequenced path of points, and
the number of subdivisions to be fit inbetween
each path-adjacent, non-extremal pair of points
(all neighboring points except the first & last),
obtain the centripetal Catmull-Rom splines that
cover each interpoint segment between neighbors.

Use these cubic polynomials to obtain coordinates
for each bounding point of an interpoint segment,
and provide them with the given points to obtain
an augmented abcissa-sequenced path of points.

The keyword `iterator`, when `false` (the default),
instructs the return of explict coordinate values.
When `iterator = true` is used, iterators over
those same coordinate values are returned.
"""
function catmullrom(points::V, nbetween::Int; 
                    augment::Bool=true, iterator::Bool=false) where {V}
    npoints = length(points)
    npoints > 3 || throw(ErrorException("four or more points are required"))
    
    if augment
        points = augment_points(points)
    end
    
    # include the given points (knots) for poly generation
    nthrough = nbetween + 2  # include both endpoints

    polys = catmullrom_polys(points)
    # interpolate using span-relative abcissae [0.0..1.0]
    xin01 = range(0.0, 1.0, length=nthrough)

    ncoords  = length(points[1])

    vals = polyval.(polys, Ref(xin01[1:end-1]))
    endval = polyval.(polys[end,:], 1.0)

    finalcoords = reshape(map(x->[x], endval), 1, ncoords)
    vals = vcat(vals, finalcoords)

    crpoints = (Iterators.flatten).(
                      [vals[:,i] for i=1:ncoords])

    return iterator ? crpoints : collect.(crpoints)
end

catmullrom(points::Array{NT,1}, nbetween::Int; augment::Bool=true, iterator::Bool=false) where {NT<:NamedTuple} =
    catmullrom(Tuple.(points), nbetween, iterator=iterator)

catmullrom(points::NTuple{N,NT}, nbetween::Int; augment::Bool=true, iterator::Bool=false) where {N, NT<:NamedTuple} =
    catmullrom(Tuple.(points), nbetween, iterator=iterator)

"""
    catmullrom_polys(points)

Traverse a sequence of points (coordinates in 2,3,..n dims),
obtaining interpolatory centripetal Catmull-Rom polynomials.
"""
function catmullrom_polys(points::V) where {V}
    npoints = length(points)
    npoints > 3 || throw(ErrorException("four or more points are required"))

    T = eltype(eltype(points))
                                 # omit the extremal points (-2)
    nspans  = npoints - 3        # count between fenceposts (-1)
    ncoords = length(points[1])  # each point has ncoordinates

    polys  = Array{Poly{T}, 2}(undef, nspans, ncoords)

    for idx = 1:nspans
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
                                pt₊::T) where {T}

    tg₀, tg₁ = centripetal_tangents(pt₋, pt₀, pt₁, pt₊)
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
                              pt₊::T) where {T}

    # intrasegment centripetal speeds
    dt₋₀ = speed(pt₋, pt₀)
    dt₀₁ = speed(pt₀, pt₁)
    dt₁₊ = speed(pt₁, pt₊)

    # correct for repeated coordinates
    ε = eps(eltype(T))
    dt₀₁ = @. ifelse(dt₀₁ < ε, one(eltype(T)), dt₀₁)
    dt₋₀ = @. ifelse(dt₋₀ < ε, dt₀₁, dt₋₀)
    dt₁₊ = @. ifelse(dt₁₊ < ε, dt₀₁, dt₁₊)

    # Compute the tangents at the two boundry points
    #   that delimit the interpolatory curve segment.
    #   - tg₀ is the tangent at pt₀
    #   - tg₁ is the tangent at pt₁

    tg₀ = @. (pt₀ - pt₋) /  dt₋₀ -
             (pt₁ - pt₋) / (dt₋₀ + dt₀₁) +
             (pt₁ - pt₀) /  dt₀₁

    tg₁ = @. (pt₁ - pt₀) /  dt₀₁ -
             (pt₊ - pt₀) / (dt₀₁ + dt₁₊) +
             (pt₊ - pt₁) /  dt₁₊

    # The tangent vectors are normalized.
    # These unit vector tangents provide
    # centripetal interpolation in [0..1].

    tg₀ = tg₀ .* dt₀₁
    tg₁ = tg₁ .* dt₀₁

    return (tg₀, tg₁)
end

"""
    centripetal_hermite(pt₀, pt₁, tg₀, tg₁)

Given the bounding points (pt₀, pt₁)
for an interpolatory curve segment
and the tangent unit vectors (tg₀, tg₁)
at those points:

Determine the coeffs (c0, c1, c2, c3)
of the cubic polynomial
   ply(t) = c0 + c1*t + c2*t^2 + c3*t^3
with associated quadratic polynomial
   dply(t) = c1 + 2*c2*t + 3*c3*t^2
such that
   ply(0) == pt₀,  dply(0) == tg₀
   ply(1) == pt₁,  dply(1) == tg₁
"""
function centripetal_hermite(pt₀::T, pt₁::T, tg₀::T, tg₁::T) where {T}
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
function speed(pt₀::T, pt₁::T) where {T}
    delta = pt₁ .- pt₀
    return sqrt(sqrtdot(delta, delta))
end

sqrtdot(a::T,b::T) where {T} = sqrt(sum(a .* b))
