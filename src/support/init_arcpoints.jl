#=
   Arcpoints:    the number of points to interpose within a delimited arc
   ArcpointsMin: the fewest    points to interpose within a delimited arc
   ArcpointsMax: the most      points to interpose within a delimited arc
=#

# make these even numbers
const defaultArcpoints = 64
const defaultArcpointsMin = 8
const defaultArcpointsMax = 256

npoints = haskey(ENV, "ARCPOINTS") ? Meta.parse(ENV["ARCPOINTS"]) : defaultArcpoints
npoints = isa(npoints, Integer) && npoints >= 0 ? npoints : defaultArcpoints
Arcpoints = npoints + isodd(npoints) # force even

npoints = haskey(ENV, "ARCPOINTS_MIN") ? Meta.parse(ENV["ARCPOINTS_MIN"]) : defaultArcpointsMin
npoints = isa(npoints, Integer) && npoints >= 0 ? npoints : defaultArcpointsMin
ArcpointsMin = npoints + isodd(npoints) # force even

npoints = haskey(ENV, "ARCPOINTS_MAX") ? Meta.parse(ENV["ARCPOINTS_MAX"]) : defaultArcpointsMax
npoints = isa(npoints, Integer) && npoints > 3 ? npoints : defaultArcpointsMax
ArcpointsMax = npoints + isodd(npoints) # force even
