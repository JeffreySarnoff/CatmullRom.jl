#  CentripetalCatmullRom.jl

### Centripetal parameterization for Catmull-Rom cubic interpolants. 


----

#### Copyright ©&thinsp;2018 by Jeffrey Sarnoff. &nbsp;&nbsp;  This work is released under The MIT License.


-----


[![][pkg-0.7-img]][pkg-0.7-url]  [![][travis-img]][travis-url]


-----

### Exports

#### catmullrom(points, interpolants)
- `points` is a tuple of points-as-tuples
- `interpolants` is a tuple of values from 0.0 to 1.0 (inclusive)
-  yields interpolating points from points[2] through points[end-1] (inclusive)


#### into01((xs...,)), into01([xs...,])
- maps values into 0.0:1.0, linearly
- minimum(xs) --> 0.0, maximum(xs) --> 1.0

```julia
julia> points2D = ([(sin(x*pi),cos(x*pi)) for x=0.0:(0.5/3):0.5]...,)
((0.0, 1.0), (0.49999999999999994, 0.8660254037844387), (0.8660254037844386, 0.5000000000000001), (1.0, 6.123233995736766e-17))

julia> polys = catmullrom_polys(points2D)
2-element Array{Poly,1}:
 Poly(0.49999999999999994 + 0.43301270189221924*x - 0.017949192431122307*x^2 - 0.049038105676658006*x^3)
 Poly(0.8660254037844387 - 0.24999999999999994*x - 0.16506350946109644*x^2 + 0.049038105676658006*x^3)  
```
-----

### Notes

With Centripetal Catmull Rom interpolation, the distances are not uniform.
Each interval is square root of the Euclidean distance between the points.

### Refs

[Parameterization and Applications of Catmull-Rom Curves](http://www.cemyuksel.com/research/catmullrom_param/catmullrom_cad.pdf)
by Cem Yuksel, Scott Schaefer, John Keyser

[The Centripetal Catmull-Rom Spline](https://howlingpixel.com/wiki/Centripetal_Catmull%E2%80%93Rom_spline)

[Catmull-Rom spline without cusps or self-intersections](https://stackoverflow.com/questions/9489736/catmull-rom-curve-with-no-cusps-and-no-self-intersections/23980479#23980479)


----

[travis-img]: https://travis-ci.org/JeffreySarnoff/CentripetalCatmullRom.jl.svg?branch=master
[travis-url]: https://travis-ci.org/JeffreySarnoff/CentripetalCatmullRom.jl



[pkg-0.6-img]: http://pkg.julialang.org/badges/CentripetalCatmullRom_0.6.svg
[pkg-0.6-url]: http://pkg.julialang.org/?pkg=CentripetalCatmullRom&ver=0.6
[pkg-0.7-img]: http://pkg.julialang.org/badges/CentripetalCatmullRom_0.7.svg
[pkg-0.7-url]: http://pkg.julialang.org/?pkg=CentripetalCatmullRom&ver=0.7
