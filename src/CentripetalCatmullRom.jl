__precompile__()

module CentripetalCatmullRom

export catmullrom,
       uniformsep, chebroots01, zero_chebroots_one, 
       into01, clamp01

using Polynomials
import Polynomials: polyval, polyder

using LinearAlgebra: dot

include("catmullrom.jl")
include("interpolant.jl")


end # module CentripetalCatmullRom
