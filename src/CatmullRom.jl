"""
    CatmullRom

Interpolate using the centripetal parameterization of Catmull-Rom splines.
This is the form that does not crimp, provides taughtness, and looks good.
"""
module CatmullRom

export catmullrom     # populates with points placed between those given

using Polynomials

# available as `CatmullRom.extendbounds`, `CatmullRom.thiele3`, etc
# reflect, linear, quadratic, thiele3, thiele4  # also for interpolation

const Type1 = AbstractArray{AbstractArray{T,1},1} where {T}
const Type2 = AbstractArray{NTuple{N,T},1} where {N,T}
const Type3 = NTuple{N, NTuple{M,T}} where {N,M,T}
const POINTS = Union{Type1, Type2, Type3} where {N,M,T}

include("centripetal/catmullrom.jl")
include("centripetal/support.jl")
include("centripetal/arcbased.jl")

include("fewpoints/twopoints.jl")
include("fewpoints/threepoints.jl")
include("fewpoints/fourpoints.jl")

end # CatmullRom
