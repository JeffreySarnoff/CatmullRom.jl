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
    return centripetal_hermite_cubic(pt₀, pt₁, tg₀, tg₁)
end

centripetal_catmullrom(pt₋::T, pt₀::T, pt₁::T, pt₊::T) where {T<:NamedTuple} =
    centripetal_catmullrom(Tuple(pt₋), Tuple(pt₀), Tuple(pt₁), Tuple(pt₊))

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
    threshold = eps(eltype(T))
    unity = one(eltype(T))
 
    dt₀₁ = @. ifelse(dt₀₁ < threshold, unity, dt₀₁)
    dt₋₀ = @. ifelse(dt₋₀ < threshold, unity, dt₀₁)
    dt₁₊ = @. ifelse(dt₁₊ < threshold, unity, dt₀₁)

    # Compute the tangents at the two boundry points
    # that delimit the interpolatory curve segment.
    #   tg₀ is the tangent at pt₀
    #   tg₁ is the tangent at pt₁

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
    centripetal_hermite_cubic(pt₀, pt₁, tg₀, tg₁)

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
function centripetal_hermite_cubic(pt₀::T, pt₁::T, tg₀::T, tg₁::T) where {T}
    c0 = pt₀
    c1 = tg₀
    c2 = @. -3*pt₀ + 3*pt₁ - 2*tg₀ - tg₁
    c3 = @.  2*pt₀ - 2*pt₁ +   tg₀ + tg₁

    return (c0, c1, c2, c3)
end

function speed(pt₀::T, pt₁::T) where {T}
    delta = pt₁ .- pt₀
    return sqrt(sqrt(dot(delta, delta)))
end
