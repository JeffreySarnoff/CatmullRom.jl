__precompile__()

module CatmullRom

export catmullrom,    # points, interpolants --> points, interpolated points
       Thiele3, Omit,
       uniformspacing, into01,
       Poly, polyval, polyder, polyint  # reexported

using Polynomials
import Polynomials: Poly, polyval, polyder, polyint

import LinearAlgebra: dot, norm

const Thiele3 = :Thiele3
const Omit    = :Omit

include("unions.jl")
include("catmullrom.jl")
include("interpolant.jl")
include("arcbased.jl")

end # module CentripetalCatmullRom
