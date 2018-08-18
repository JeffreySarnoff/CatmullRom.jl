"""
   CatmullRom

This package provides open and closed centripetal Catmull-Rom splines.
The centripetal parameterization of Catmull-Rom splines works best in use:
it does not crimp and provides a useful amount of relative taughtness.

see @ref(catmullrom) for use.
"""
module CatmullRom

export catmullrom,    # points, interpolants --> points, interpolated points
       Omit, Linear, Quadratic, Thiele3,
       uniform01, into01, clamp01,
       Poly, polyval, polyder, polyint  # reexported

using Polynomials
import Polynomials: Poly, polyval, polyder, polyint

import LinearAlgebra: dot, norm

include("consts.jl")
include("unions.jl")
include("centripetal_catmullrom.jl")
include("interpolant.jl")

end # module CentripetalCatmullRom
