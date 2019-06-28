"""
    outside(points [; scale=1/16])

Determine the two extremal points (outside of the boundry points)
used to anchor the sequence of centripetal Catmull-Rom spline spans.
Incorporate these boundry points in the point sequence.
"""
function outside(points::Points; scale=1/16)
    npoints(points) > 3 || throw(DomainError("4 or more points are required"))
    firstpoints = points[1:4]
    lastpoints  = points[end-3:end]
    initialpoint, finalpoint = pointbefore(firstpoints, scale), pointafter(lastpoints, scale)
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
