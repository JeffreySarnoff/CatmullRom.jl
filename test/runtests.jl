using CentripetalCatmullRom

if VERSION >= v"0.7-"
    using Test
else
    using Base.Test
end

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

