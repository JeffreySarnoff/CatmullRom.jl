"""
    extendbounds(points [; scale=1/16])

Extend the point sequence by one point before the start and
by one point after the end. The extended point sequence
anchors the centripetal Catmull-Rom splines.
"""
function extendbounds(points::Points; scale=ReflectionScale)
    npoints(points) > 3 || throw(DomainError("4 or more points are required"))
    points = outside(points, scale=scale)
    return points
end

"""
    outside(points [; scale])

Determine the two extremal points (outside of the boundary points)
used to anchor the sequence of centripetal Catmull-Rom spline spans.
Incorporate these boundary points in the point sequence.
"""
function outside(points::Points; scale=ReflectionScale)
    initialpoint = pointbefore(points[1:4], scale)
    finalpoint   = pointafter(points[end-3:end], scale)
    pushfirst!(push!(points, finalpoint), initialpoint)
    return points
end

function pointbefore(firstpoints::Points, scale)
    initialx = reflectback(first.(firstpoints[1:2])..., scale=scale)
    initialys = thiele4(firstpoints..., initialx)
    initialpoint = (initialx, initialys...,)
    return convert(eltype(firstpoints), initialpoint)
end

function pointafter(lastpoints::Points, scale)
    finalx = reflectforward(first.(lastpoints[end-1:end])..., scale=scale)
    finalys = thiele4(lastpoints..., finalx)
    finalpoint = (finalx, finalys...,)
    return convert(eltype(lastpoints), finalpoint)
end
