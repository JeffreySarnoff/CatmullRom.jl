using CatmullRom
using Test

@testset "clamp01" begin
    for (X,Z) in ( (0.0, 0.0), (1.0f0, 1.0f0),
                   (prevfloat(1.0), prevfloat(1.0)), (nextfloat(1.0), 1.0),
                   (prevfloat(0.0), 0.0), (nextfloat(0.0), nextfloat(0.0)) 
                 )
        @eval @test clamp01($X) === $Z
    end 
end

@testset "uniform01" begin
    for (X,Z) in ( (2, [0.0, 1.0]),
                   (3, [0.0, 1.0/2.0, 1.0]),
                   (4, [0.0, 1.0/3.0, 2.0/3.0, 1.0]),
                   (5, [0.0, 1.0/4.0, 2.0/4.0, 3.0/4.0, 1.0])
                 )
        @eval uniform01($X) === $Z
    end
end

