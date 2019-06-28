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
    initialx = firstpoints[1][1] + (firstpoints[1][1] - firstpoints[2][1]) * scale
    initialys = thiele4(firstpoints..., initialx)
    initialpoint = (initialx, initialys...,)
    return convert(eltype(firstpoints), initialpoint)
end

function pointafter(lastpoints::Points, scale)
    finalx = lastpoints[end][1] + (lastpoints[end][1] - lastpoints[end-1][1]) * scale
    finalys = thiele4(lastpoints..., finalx)
    finalpoint = (finalx, finalys...,)
    return convert(eltype(lastpoints), finalpoint)
end
