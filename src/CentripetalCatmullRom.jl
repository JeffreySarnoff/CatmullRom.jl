__precompile__()

module CentripetalCatmullRom

export CatmullRom, 
       PT1D, PT2D, PT3D, PT4D, coords,   # points in 1D, 2D, 3D, 4D [coordinates :x, :y, :z, :t]
       xcoord, ycoord, zcoord, tcoord,   # access point coordinates in 1D, 2D, 3D, 4D
       Δpoint, Δpoint2, dpoint, dpoint2, # dpoint, dpoint2 alias Δpoint, Δpoint2   
       polyval, polyder                  # from Polynomials, specialized


import Base: values                      # for NamedTuples

using Polynomials
import Polynomials: polyval, polyder

include("points_1Dto4D.jl")
include("polys_1Dto4D.jl")




function test()
    p0 = POINT2D(0.0, 0.0)
    p1 = POINT2D(1.0, 1.0)
    p2 = POINT2D(1.5, 1.25)
    p3 = POINT2D(2.0, 0.0)

    crpoly = CatmullRomPolys(p0, p1, p2, p3)

    # p1 ... n-1 points ... p2
    
    n = 4
    coords = Vector{POINT2D}(undef, n+1)
    m = inv(n)

    for i=0:n
       xcoord, ycoord = polyval(crpoly, min(1.0,m*i))
       coord = POINT2D(xcoord, ycoord)
       coords[i+1] = coord
    end
    return coords
end

print(test())



end # module
