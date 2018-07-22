__precompile__()

module CentripetalCatmullRom

export catmullrom,                      # points, interpolants --> points, interpolated points
       catmullrom_onpath,               # rough and inexact pathlengths (arc traversals)
                                        #    covering curvilinear segments between points 
       uniformsep,
       chebyshevsep,
       into01,
       clamp01,
       catmullrom_polys,                 # Catmull-Rom centripetal cubics
       catmullrom_polys_d01,             #     cubics and quadratic first derivatives
       catmullrom_polys_d012,            #     cubics and quadratic first and monic second derivatives
       catmullrom_polys_i01,             #     cubics and quadrics of integration
       Poly, polyval, polyder, polyint   # reexported

using Polynomials
import Polynomials: Poly, polyval, polyder, polyint

import LinearAlgebra: dot, norm

include("catmullrom.jl")
include("interpolant.jl")
include("arcbased.jl")

end # module CentripetalCatmullRom
