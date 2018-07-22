#=
    given some adjacency-ordered ND points (at least four)
    obtain rough estimates of the path-traversal `time`
    over each curvilinear segment
    (bounded at start and at end by two [adjacent] points).
    n points implies n-1 adjacent interpoint segments.
=#



@inline function linearsep(pointa::P, pointb::P) where {N,T, P<:NTuple{N,T}}
    sqrt(lawofcosines(norm(pointa), anglesep(pointa, pointb), norm(pointb)))
end

# relatively fast determination of angular separation
#    UNCHECKED PRECONDITION:
#       both points are given relative to the same origin
#
#  >>>  for a numerically rigourous approach to angular separation
#  >>>  use AngleBetweenVectors.jl
#

function anglesep(pointa::NTuple{N,T}, pointb::NTuple{N,T}) where {N,T}
    dota = dot(pointa, pointa)
    dotb = dot(pointb, pointb)
    iszero(dota) || iszero(dotb) && return zero(T)
    
    dotb = sqrt(dota * dotb)
    dota = dot(pointa, pointb)
    acos( dota / dotb )
end
    
@inline lawofcosines(side1, angle2sides, side2) =
    side1*side1 + side2*side2 - side1*side2 * 2*cos(angle2sides)
   

    
#=
    Given 4 ND points, roughly approximate the arclength
    of the centripetal Catmull-Rom curvilinear segment 
    that would be determined by two bounding points
    and the tangents they determine.
    
    this algorithm was developed by Jens Gravesen
=#
function rough_trisegment_arclength(pts::NTuple{4, NTuple{N,T}}) where {N, T}
     ldist12 = linearsep(pts[2], pts[1])
     ldist23 = linearsep(pts[3], pts[2])
     ldist34 = linearsep(pts[4], pts[3])
     ldist14 = linearsep(pts[4], pts[1])
     
     linesegments = ldist12 + ldist23 + ldist34
     arclength = (linesegments + ldist14) / 2
     # estimated_error  = linesegments - ldist14
     
     return arclength   
end

#=
    rough approximation to the length of the arc
       connecting points 2 and 3, this is center arc
       that is interpolated between with each use of
       catmullrom_4points 
=#

function rough_midsegment_arclength(pts::NTuple{4, NTuple{N,T}}) where {N, T}
     ldist12 = linearsep(pts[2], pts[1])
     ldist23 = linearsep(pts[3], pts[2])
     ldist34 = linearsep(pts[4], pts[3])
     ldist14 = linearsep(pts[4], pts[1])
     
     linesegments = ldist12 + ldist23 + ldist34
     midsegment_proportionalweight  = ldist23 / linesegments

     arclength = (linesegments + ldist14) / 2
     midsegment_arclength  = arclength * midsegment_proportionalweight
    
     return midsegment_arclength
end

