using CatmullRom
using Test

@testset "closed curve interpositions" begin
  xs = [sinpi(t) for t=range(0.0, stop=2.0, length=5)];
  ys = [cospi(t) for t=range(0.0, stop=2.0, length=5)];
  points = collect(zip(xs,ys));
  cr_xys = catmullrom(points, 2)
  @test length(cr_xys[1]) == length(xs) + 2*(length(xs)-1)
  @test length(points) == length(xs) + 2
end

#=
    example from
    Block-based Thiele-like blending rational interpolation
    Qian-Jin Zhao, Jieqing Tan
    Journal of Computational and Applied Mathematics
    https://doi.org/10.1016/j.cam.2005.03.089
=#
@testset "interpolation" begin
  xs=[0.0, 1.0, 2.0,  3.0,  4.0, 5.0]
  ys=[0.0, 1.0, 4.0, -1.0, -2.0, 2.0]
  points = collect(zip(xs,ys))
  @test CatmullRom.thiele4(points[1:4]..., 6.0)[1] == 2.0 + 2/3
  @test CatmullRom.thiele4(points[2:5]..., 1.0)[1] == 1.0
  @test CatmullRom.quadratic(points[3:5]..., 1.5)[1] == 8.0
end

@testset "approximate arclengths"
  xs=[0.0, 1.0, 2.0,  3.0,  4.0, 5.0]
  ys=[0.0, 1.0, 4.0, -1.0, -2.0, 2.0]
  points = collect(zip(xs,ys))
  @test length(points) - 2*2 + 1 == length(CatmullRom.approx_catmullrom_arclengths(points))  
end
