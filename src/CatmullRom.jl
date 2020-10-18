"""
    CatmullRom

Interpolate using the centripetal parameterization of Catmull-Rom splines.
This is the form that does not crimp, provides taughtness, and looks good.

See: [`catmullrom`](@ref), [`catmullrom_by_arclength`](@ref)

See also: [`close_seq!`](@ref)
""" CatmullRom

module CatmullRom

export catmullrom,               # populates with points placed between those given
       catmullrom_by_arclength,  # like catmullrom with points placed relative to arclength
       close_seq!,               # ensures a point sequence is closed, precisely
       floatvecs                 # obtain a vector of vectors of floats

using LinearAlgebra: dot, norm, normalize
using Polynomials: Poly, polyval

using LoopVectorization, StructArrays, StaticArrays

# include("support/basemethods.jl")
include("support/init_arcpoints.jl")

include("point/coordinate.jl")
include("point/sequence.jl")

include("centripetal/approx_arclength.jl")

# The sorts of sequences understood to hold point coordinates
#     defines `Points`m `npoints(Points)`, `ncoords(Points)`
include("point/sequences.jl")

# suggest outermost two points for Catmull-Rom spline sequence
const ReflectionScale = 1.0

include("fewpoints/outside.jl")

# centripetal Catmull-Rom compution
include("centripetal/catmullrom1.jl")
include("centripetal/catmullrom2.jl")
include("centripetal/arcbased.jl")

# extrapolation using low-order interpolatory fits
include("fewpoints/twopoints.jl")    # linear fit, reflection
include("fewpoints/threepoints.jl")  # thiele3 rational fit, general quadratic fit
include("fewpoints/fourpoints.jl")   # thiele4 rational fit (preferred)

# myunzip, anglesep
include("presence/support.jl")

end # CatmullRom
