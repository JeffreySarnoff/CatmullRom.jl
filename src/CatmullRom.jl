"""
    CatmullRom

Interpolate using the centripetal parameterization of Catmull-Rom splines.
This is the form that does not crimp, provides taughtness, and looks good.
"""
module CatmullRom

export catmullrom     # populates with points placed between those given

using LinearAlgebra: dot, norm, normalize

using Polynomials: Poly, polyval, polyder, polyint


# The sorts of sequences understood to hold point coordinates
#     defines `Points`m `npoints(Points)`, `ndims(Points)`
include("pointsequences.jl")


# centripetal Catmull-Rom compution
include("centripetal/catmullrom.jl")
include("centripetal/support.jl")
include("centripetal/arcbased.jl")


# suggest outermost two points for Catmull-Rom spline sequence
include("fewpoints/outside.jl")

# extrapolation using low-order interpolatory fitting
include("fewpoints/twopoints.jl")    # linear fit, reflection
include("fewpoints/threepoints.jl")  # thiele3 rational fit, general quadratic fit
include("fewpoints/fourpoints.jl")   # thiele4 rational fit (preferred)


end # CatmullRom
