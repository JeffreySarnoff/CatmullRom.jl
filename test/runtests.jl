using CentripetalCatmullRom
using Test

@test true

#=
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
=#

