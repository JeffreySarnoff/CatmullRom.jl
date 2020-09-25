const defaultArcpoints = 64
const defaultArcpointsMin = 8
const defaultArcpointsMax = 256

npoints = haskey(ENV, "ARCPOINTS") ? Meta.parse(ENV["ARCPOINTS"]) : defaultArcpoints
const Arcpoints = isa(npoints, Integer) ? defaultArcpoints : npoints

npoints = haskey(ENV, "ARCPOINTS_MIN") ? Meta.parse(ENV["ARCPOINTS_MIN"]) : defaultArcpointsMin
const ArcpointsMin = isa(npoints, Integer) ? defaultArcpointsMin : npoints

npoints = haskey(ENV, "ARCPOINTS_MAX") ? Meta.parse(ENV["ARCPOINTS_MAX"]) : defaultArcpointsMax
const ArcpointsMax = !isa(npoints, Integer) ? defaultArcpointsMax : npoints
