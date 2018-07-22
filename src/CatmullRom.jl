__precompile__()

module CatmullRom

export catmullrom,                      # points, interpolants --> points, interpolated points

       catmullrom_extents,              # rough and inexact pathlengths (arc traversals)
                                        #    covering curvilinear segments between points 
       uniformsep,                      # gen your interpolant values for uniform spacing 
       chebyshevsep,                    # gen your interpolant values for chebyshev spacing

       Poly, polyval, polyder, polyint  # reexported

using Polynomials
import Polynomials: Poly, polyval, polyder, polyint

import LinearAlgebra: dot, norm


# D is the dimensionality of the coordinate space
# R is the eltype of points, the numeric type underlying each coordinate
# M is the number of points, when using a tuple of points
const TupledPoint{D,R} = NTuple{D,R} where {D,R}  
const TupledPoints{M,D,R} = NTuple{M, TupledPoint{D,R}} where {M,D,R}
# F is the numeric type used for values (usually should match R above)
# L is the number of values, when using a tuple
const ValueSeq{L,T} = Union{AbstractArray{T,1}, NTuple{L,T}} where {L,T}
const PointSeq{M,D,T} = Union{AbstractArray{T,1}, TupledPoints{M,D,T}} where {M,D,T}

include("catmullrom.jl")
include("interpolant.jl")
include("arcbased.jl")

end # module CentripetalCatmullRom
