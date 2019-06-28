"""
    CatmullRom

Interpolate using the centripetal parameterization of Catmull-Rom splines.
This is the form that does not crimp, provides taughtness, and looks good.
"""
module CatmullRom

export catmullrom     # populates with points placed between those given

#=
    normalize(v) so that its p-norm equals unity, i.e. norm(v, p) == 1. 
=#
using LinearAlgebra: dot, norm, normalize

#=
    Poly([coeff_0, coeff_1, .., coeff_i, .., coeff_n])
          coeff_0 is the constant term
          coeff_i is the coefficient of x^i
          coeff_n is the coefficient of x^n for a polynomial of degree n

   polyval(apoly, atvalue)
   polyder(apoly, nth)      nth derivative of apoly
   polyint(apoly, k=0)      integral of apoly with arbitrary constant k added
=#
using Polynomials: Poly, polyval, polyder, polyint


# available as `CatmullRom.extendbounds`, `CatmullRom.thiele3`, etc
# reflect, linear, quadratic, thiele3, thiele4  # also for interpolation

# The sorts of sequences understood to hold point coordinates
#     defines `Points`m `npoints(Points)`, `ndims(Points)`
include("pointsequences.jl")

include("centripetal/catmullrom.jl")
include("centripetal/support.jl")
include("centripetal/arcbased.jl")

include("fewpoints/twopoints.jl")
include("fewpoints/threepoints.jl")
include("fewpoints/fourpoints.jl")
include("fewpoints/outside.jl")

end # CatmullRom
