__precompile__()

module CentripetalCatmullRom

export CatmullRom, 
       PT1D, PT2D, PT3D, PT4D, coords,   # points in 1D, 2D, 3D, 4D [coordinates :x, :y, :z, :t]
       xcoord, ycoord, zcoord, tcoord,   # access point coordinates in 1D, 2D, 3D, 4D
       Δpoint, Δpoint2, dpoint, dpoint2, # dpoint, dpoint2 alias Δpoint, Δpoint2   
       polyval, polyder,                 # from Polynomials, specialized
       δpolyval, dpolyval                # dpolyval aliases δpolyval 

import Base: values                      # for NamedTuples

using Polynomials
import Polynomials: polyval, polyder

include("points_1Dto4D.jl")
include("polys_1Dto4D.jl")

end # module CentripetalCatmullRom
