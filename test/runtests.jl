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

