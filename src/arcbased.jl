function catmullrom_pathparts(points::Vector{NTuple{N,T}}; subdivisions::Int=64) where {N,T}
    nspans = length(points) - 1
    # each span gets at least 3: 0.0, mid, 1.0; -1 corrects overlap
    ninterpolants  = (nspans * (subdivisions + 2)) - 1
    spancounts = catmullrom_onpath(points, ninterpolants)
    return spancounts
end    

function catmullrom_onpath(points::Vector{NTuple{N,T}}, ninterpolants::Int) where {N,T}
    extents = catmullrom_extents(points)
    # use extents to apportion ninterpolants
    relspans = extents .* inv(sum(extents))    # relspans sum to 1
    spancounts = trunc.(Int, round.((relspans .* ninterpolants), RoundNearest))
    
    while sum(spancounts) > ninterpolants
        idx = findfirst(spancounts .== maximum(spancounts))
        spancounts[idx] -= 1
    end
    while sum(spancounts) < ninterpolants
        idx = findfirst(spancounts .== minimum(spancounts))
        spancounts[idx] += 1
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

function catmullrom_extents(points::Vector{NTuple{N,T}}) where {N,T}
    npoints = length(points)
    result = Vector{T}(npoints - 1)
    result[1]   = linearsep(points[1], points[2])
    result[end] = linearsep(points[npoints-1], points[npoints])
    
    ngroupsof4 = npoints - 3
    for idx in 1:ngroupsof4
        fourpoints = (points[idx:idx+3]...,)
        arclength = rough_centralsegment_arclength(fourpoints)
        result[idx+1] = arclength
    end

    return result
end


        
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
    (iszero(dota) || iszero(dotb)) && return zero(T)
    
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

function rough_centralsegment_arclength(pts::NTuple{4, NTuple{N,T}}) where {N, T}
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

#=
0.02245448808080272	5.767944580792328e-6	1.0005138773654616
0.11034730555433708	3.492803405297229e-5	1.000633256844574
0.20680347495224183	0.0007521441574384058	1.0073005513192965
0.27722816617416257	0.0007397400132340692	1.0053509654889026
0.02953427762165981	0.0003211994629658821	1.0219901142372627
0.02836586860799523	0.0010214234384144927	1.0747079293128807
0.22249019182783417	0.031536873165516120	1.3303097677111957
0.28190746915694587	0.030669273789136087	1.244144993512922
0.01749704778645134	0.007188486833679952	2.3946634969932195
0.1082295845425749	5.318705844575833e-5	1.0009833394286136

julia> function rough_centralsegment_arclength(pts::NTuple{4, NTuple{N,T}}) where {N, T}
            ldist12 = linearsep(pts[2], pts[1])
            ldist23 = linearsep(pts[3], pts[2])
            ldist34 = linearsep(pts[4], pts[3])
            ldist14 = linearsep(pts[4], pts[1])
            
            linesegments = ldist12 + ldist23 + ldist34
            midsegment_proportionalweight  = ldist23 / linesegments

            arclength = (linesegments + ldist14) / 2
            midsegment_arclength  = arclength * midsegment_proportionalweight
         estimated_error  = linesegments - ldist14  ; println(midsegment_arclength,"\t", estimated_error, "\t",inv((midsegment_arclength-estimated_error)/(midsegment_arclength+estimated_error)))
            return midsegment_arclength
       end
=#
