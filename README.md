# CatmullRom.jl

### Centripetal parameterization for Catmull-Rom cubic interpolants. 


#### Copyright ©&thinsp;2018-2019 by Jeffrey Sarnoff. &nbsp;&nbsp;  This work is released under The MIT License.


-----


[![Build Status](https://travis-ci.org/JeffreySarnoff/CatmullRom.jl.svg?branch=master)](https://travis-ci.org/JeffreySarnoff/CatmullRom.jl)&nbsp;&nbsp;&nbsp;[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](http://jeffreysarnoff.github.io/CatmullRom.jl/stable/)&nbsp;&nbsp;&nbsp;[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](http://jeffreysarnoff.github.io/CatmullRom.jl/dev/)&nbsp;&nbsp;&nbsp;[![codecov](https://codecov.io/gh/JeffreySarnoff/CatmullRom.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JeffreySarnoff/CatmullRom.jl)


-----

### Exports

## constructor

- catmullrom(points, how_many_interpolants)
    - `points` is a vector or tuple of points-as-tuples or points-as-vectors
    - `how_many_interpolants` counts the points to be interpolated
    -  yields interpolating points from points[2] through points[end-1] (inclusive)

### closed curves

To close a curve make sure that the first point and the last point are the same point.


```julia

julia> points = ([(sinpi(x),cospi(x)) for x=0.0f0:(0.25f0/3.0f0):0.25f0]...,)
((0.0f0, 1.0f0), (0.25881904f0, 0.9659258f0), (0.5f0, 0.8660254f0), (0.70710677f0, 0.70710677f0))

julia> xs,ys = CatmullRom.unzip(points)
4×2 Array{Float32,2}:
 0.0       1.0
 0.258819  0.965926
 0.5       0.866025
 0.707107  0.707107

julia> crpoints = catmullrom(points, 2)
2-element Array{Array{Float32,1},1}:
 [0.25881904, 0.34178123, 0.4227836, 0.5]
 [0.9659258, 0.93968755, 0.9061352, 0.8660252]


 
```
-----

### Notes

With Centripetal Catmull Rom interpolation, the distances are not uniform.
Each interval is square root of the Euclidean distance between the points.

----

### Refs

[Parameterization and Applications of Catmull-Rom Curves](http://www.cemyuksel.com/research/catmullrom_param/catmullrom_cad.pdf)

[The Centripetal Catmull-Rom Spline](https://howlingpixel.com/wiki/Centripetal_Catmull%E2%80%93Rom_spline)

[Catmull-Rom spline without cusps or self-intersections](https://stackoverflow.com/questions/9489736/catmull-rom-curve-with-no-cusps-and-no-self-intersections/23980479#23980479)

-----

[travis-img]: https://travis-ci.org/JeffreySarnoff/CatmullRom.jl.svg?branch=master
[travis-url]: https://travis-ci.org/JeffreySarnoff/CatmullRom.jl

[pkg-1.0-img]: http://pkg.julialang.org/badges/CatmullRom_1.0.svg
[pkg-1.0-url]: http://pkg.julialang.org/?pkg=CatmullRom&ver=1.0
