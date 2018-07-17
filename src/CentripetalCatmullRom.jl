__precompile__()

module CentripetalCatmullRom

export catmullrom,
       uniformsep, chebroots01, zero_chebroots_one, 
       into01, clamp01,
       Poly, polyval, polyder, polyint

using Polynomials
import Polynomials: Poly, polyval, polyder, polyint

import LinearAlgebra: dot

include("catmullrom.jl")
include("centripetal.jl")
include("interpolant.jl")

end # module CentripetalCatmullRom
