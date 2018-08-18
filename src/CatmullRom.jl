module CatmullRom

export catmullrom,    # points, interpolants --> points, interpolated points
       Omit, Linear, Quadratic, Thiele3,
       uniformspacing, into01, clamp01,
       Poly, polyval, polyder, polyint  # reexported

using Polynomials
import Polynomials: Poly, polyval, polyder, polyint

import LinearAlgebra: dot, norm

include("consts.jl")
include("unions.jl")
include("centripetal_catmullrom.jl")
include("interpolant.jl")

end # module CentripetalCatmullRom
