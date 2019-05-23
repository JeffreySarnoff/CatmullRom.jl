function catmullrom_pathparts(points::P, avg_points_per_segment::Int=19) where {P}
    points_to_interpolate = avg_points_per_segment * (length(points) - 1)
    spancounts = catmullrom_onpath(points, points_to_interpolate)
    return spancounts
end

function catmullrom_onpath(points::P, ninterpolants::Int) where {P}
    extents = catmullrom_extents(points)
    # use extents to apportion ninterpolants
    relspans = extents .* inv(sum(extents))    # relspans sum to 1
    spancounts = trunc.(Int, round.((relspans .* ninterpolants), RoundNearest))

    if any(iszero.(spancounts))
        spancounts = spancounts .+ 1
    end

    while sum(spancounts) < ninterpolants
        idx = findfirst(spancounts .== minimum(spancounts))
        spancounts[idx] += 1
        sum(spancounts) < ninterpolants &&
           let idx = findlast(spancounts .== minimum(spancounts));
               spancounts[idx] += 1
           end
    end
    while sum(spancounts) > ninterpolants
        idx = findfirst(spancounts .== maximum(spancounts))
        spancounts[idx] -= 1
        sum(spancounts) > ninterpolants &&
           let idx = findlast(spancounts .== maximum(spancounts));
               spancounts[idx] -= 1
           end
    end

    return spancounts
end


#=
    given some adjacency-ordered ND points (at least four)
    obtain rough estimates of the path-traversal `time`
    over each curvilinear segment
    (bounded at start and at end by two [adjacent] points).
    n points implies n-1 adjacent interpoint segments.
=#

function catmullrom_extents(points::P) where {p}
    npoints = length(points)
    result = Vector{R}(undef, npoints - 1)
    result[1]   = linearseparation(points[1], points[2])
    result[end] = linearseparation(points[npoints-1], points[npoints])

    ngroupsof4 = npoints - 3
    for idx in 1:ngroupsof4
        fourpoints = (points[idx:idx+3]...,)
        arclength = approximate_arclength(fourpoints) #  for central segment
        result[idx+1] = arclength
    end

    return result
end


#=
    Given 4 ND points, roughly approximate the arclength
    of the centripetal Catmull-Rom curvilinear segment
    that would be determined by two bounding points
    and the tangents they determine.
    this algorithm was developed by Jens Gravesen
=#
function approximate_arclength(points::P) where {P}
     ldist12 = linearseparation(points[2], points[1])
     ldist23 = linearseparation(points[3], points[2])
     ldist34 = linearseparation(points[4], points[3])
     ldist14 = linearseparation(points[4], points[1])

     linesegments = ldist12 + ldist23 + ldist34
     arclength = (linesegments + ldist14) * ldist23
     arclength /= 2 * linesegments

     # arclength = (linesegments + ldist14) / 2
     # arclength *=  ldist23 / linesegments
     # errorest  = linesegments - ldist14

     return arclength
end


linearseparation(a::P1, b::P1) where {P1} =
    sqrt(lawofcosines(norm(a), anglesep(a, b), norm(b)))

lawofcosines(side1, anglebetween, side2) =
    side1*side1 + side2*side2 - side1*side2 * 2*cos(anglebetween)


# relatively fast determination of angular separation
#    UNCHECKED PRECONDITION:
#       both points are given relative to the same origin
#
#  >>>  for a numerically rigourous approach to angular separation
#  >>>  use AngleBetweenVectors.jl
#

function anglesep(pointa::P1, pointb::P1) where {P1}
    dota = dot(pointa, pointa)
    dotb = dot(pointb, pointb)
    (iszero(dota) || iszero(dotb)) && return zero(T)

    dotb = sqrt(dota * dotb)
    dota = dot(pointa, pointb)
    acos( dota / dotb )
end


#=
    rough approximation to the length of the arc
       connecting points 2 and 3, this is center arc
       that is interpolated between with each use of
       catmullrom_4points
=#
#=
function rough_centralsegment_arclength(points::PointSeq{M,D,R}) where {M,D,R}
     ldist12 = linearsep(points[2], points[1])
     ldist23 = linearsep(points[3], points[2])
     ldist34 = linearsep(points[4], points[3])
     ldist14 = linearsep(points[4], points[1])

     linesegments = ldist12 + ldist23 + ldist34
     midsegment_proportionalweight  = ldist23 / linesegments
     arclength = (linesegments + ldist14) / 2
     midsegment_arclength  = arclength * midsegment_proportionalweight

     return midsegment_arclength
end
=#
