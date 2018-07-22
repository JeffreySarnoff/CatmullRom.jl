__precompile__()

module CatmullRom

export catmullrom    # points, interpolants --> points, interpolated points
#=
       catmullrom_extents,              # rough and inexact pathlengths (arc traversals)
                                        #    covering curvilinear segments between points 
       uniformsep,                      # gen your interpolant values for uniform spacing 
       chebyshevsep,                    # gen your interpolant values for chebyshev spacing

       Poly, polyval, polyder, polyint  # reexported
=#

using Polynomials
import Polynomials: Poly, polyval, polyder, polyint

import LinearAlgebra: dot, norm


# D is the dimensionality of the coordinate space
# R is the eltype of points, the numeric type underlying each coordinate
# M is the number of points, when using a tuple of points
# F is the numeric type used for values (usually should match R above)
# L is the number of values, when using a tuple


const TupAsPoint{D,R} = NTuple{D,R} where {D,R}
const VecAsPoint{R}   = Vector{R}   where {R}
const OnePoint{D,R}   = Union{TupAsPoint{D,R}, VecAsPoint{R}}  where {D,R}

const TupAsTupPoints{M,D,R} = NTuple{M, TupAsPoint{D,R}} where {M,D,R}
const TupAsVecPoints{M,R}   = NTuple{M, VecAsPoint{R}}   where {M,R}
const VecAsTupPoints{D,R}   = Vector{TupAsPoint{D,R}}    where {D,R}
const VecAsVecPoints{R}     = Vector{VecAsPoint{R}}      where {R}

const PointSeq{M,D,R}  = Union{TupAsTupPoints{M,D,R}, TupAsVecPoints{M,R}, 
                               VecAsTupPoints{D,R}, VecAsVecPoints{R}} where {M,D,R}

const TupAsValues{L,F} = NTuple{L,F} where {L,F}
const VecAsValues{F}   = Vector{F} where {F}
const ValueSeq{L,F}    = Union{TupAsValues{L,F}, VecAsValues{F}} where {L,F}


include("catmullrom.jl")
include("interpolant.jl")
include("arcbased.jl")

end # module CentripetalCatmullRom
