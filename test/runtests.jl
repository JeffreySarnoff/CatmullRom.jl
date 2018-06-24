using CentripetalCatmullRom
using LinearAlgebra
using Test

points2D = ((0.0, 0.0), (0.25, 2.0), (0.75, 2.0), (1.0, 0.0))
points3D = ((0.0, 0.0, 0.0), (0.25, 2.0, 0.5), (1.75, 2.0, 1.0), (1.0, 0.0, 0.0))

interpolants = collect(0.0:0.2:1.0)

pts2D = catmullrom(points2D, interpolants)
pts3D = catmullrom(points3D, interpolants)

@test Float32.(pts2D[1,:]) == Float32.([0.25, 2.0])
@test Float32.(pts3D[1,:]) == Float32.([1.75, 0.25, 2.0])


@test true

#=
function test()
    p0 = PT2D(0.0, 0.0)
    p1 = PT2D(1.0, 1.0)
    p2 = PT2D(1.5, 1.25)
    p3 = PT2D(2.0, 0.0)

    crcubics = CatmullRomCubics(p0, p1, p2, p3)

    # p1 ... n-1 points ... p2
    
    n = 4
    coords = Vector{PT2D}(undef, n+1)
    m = inv(n)

    for i=0:n
       xcoord, ycoord = polyval(crcubics, min(1.0,m*i))
       coord = PT2D(xcoord, ycoord)
       coords[i+1] = coord
    end
    return coords
end
=#

