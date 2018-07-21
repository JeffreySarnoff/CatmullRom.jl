__precompile__()

module CentripetalCatmullRom

export catmullrom,
       uniformsep, chebyshevsep, 
       into01, clamp01,
       Poly, polyval, polyder, polyint

using Polynomials
import Polynomials: Poly, polyval, polyder, polyint

import LinearAlgebra.dot

include("catmullrom.jl")
include("interpolant.jl")

end # module CentripetalCatmullRom
