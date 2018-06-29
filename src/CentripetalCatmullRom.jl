__precompile__()

module CentripetalCatmullRom

export catmullrom,
       uniformsep, chebroots01, zero_chebroots_one, 
       into01, clamp01,
       Poly, polyval, polyder, polyint

using Polynomials
import Polynomials: Poly, polyval, polyder, polyint

using LinearAlgebra: dot

include("catmullrom.jl")
include("cetripetal.jl")
include("interpolant.jl")

end # module CentripetalCatmullRom
