__precompile__()

module CatmullRom

export catmullrom,    # points, interpolants --> points, interpolated points
       Omit, Linear, Quadratic, Thiele3,
       uniformspacing, into01,
       Poly, polyval, polyder, polyint  # reexported

using Polynomials
import Polynomials: Poly, polyval, polyder, polyint

import LinearAlgebra: dot, norm

const Quadratic = :Quadratic
const Thiele3   = :Thiele3
const Omit      = :Omit

include("unions.jl")
include("catmullrom.jl")
include("interpolant.jl")

end # module CentripetalCatmullRom
