__precompile__()

module CentripetalCatmullRom

export catmullrom,
       uniformsep, chebyshevsep,
       into01, clamp01,
       ccr_polys, ccr_polys_dpolys, ccr_polys_ipolys, ccr_allpolys,
       Poly, polyval, polyder, polyint

using Polynomials
import Polynomials: Poly, polyval, polyder, polyint

import LinearAlgebra: dot, norm

include("catmullrom.jl")
include("interpolant.jl")

end # module CentripetalCatmullRom
