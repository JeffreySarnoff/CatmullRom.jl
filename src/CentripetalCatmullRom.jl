__precompile__()

module CentripetalCatmullRom

export catmullrom,                      # points, interpolants --> points, interpolated points

       catmullrom_extents,              # rough and inexact pathlengths (arc traversals)
                                        #    covering curvilinear segments between points 
       uniformsep,                      # gen your interpolant values for uniform spacing 
       chebyshevsep,                    # gen your interpolant values for chebyshev spacing

       Poly, polyval, polyder, polyint  # reexported

using Polynomials
import Polynomials: Poly, polyval, polyder, polyint

import LinearAlgebra: dot, norm



include("catmullrom.jl")
include("interpolant.jl")
include("arcbased.jl")

end # module CentripetalCatmullRom
