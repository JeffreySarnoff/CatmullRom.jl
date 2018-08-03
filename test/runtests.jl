using CatmullRom
using Test

points2D = ((0.0, 0.0), (0.25, 2.0), (0.75, 2.0), (1.0, 0.0))
points3D = ((0.0, 0.0, 0.0), (0.25, 2.0, 0.5), (1.75, 2.0, 1.0), (1.0, 0.0, 0.0))

interpolants = collect(0.0:0.2:1.0)

pts2D = catmullrom(points2D, interpolants, endpoints=Omit)
pts3D = catmullrom(points3D, interpolants, endpoints=Omit)

@test Float32.(pts2D[1,:][1]) == Float32.((0.25, 2.0))
@test Float32.(pts3D[1,:][1]) == Float32.((0.25, 2.0, 0.5))


points2D = ((0.0, 0.0), (0.25, 2.0), (0.5, 3.0), (0.75, 2.0), (1.0, 0.0))
interpolants = collect(0.0:0.2:1.0)
pts2D = catmullrom(points2D, interpolants, endpoints=Omit)
@test pts2D[6,:][1] == (0.5, 3.0)
@test Float32.(pts2D[5,:][1]) == (0.45015842f0, 2.9427261f0)
