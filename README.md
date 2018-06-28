#  CentripetalCatmullRom.jl
Centripetal variant


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

-----

### Notes

With Centripetal Catmull Rom interpolation, the distances are not uniform.
Each interval is square root of the Euclidean distance between the points.

### Refs

Parameterization and Applications of Catmull-Rom Curves
Cem Yuksel, Scott Schaefer, John Keyser
http://www.cemyuksel.com/research/catmullrom_param/catmullrom_cad.pdf


https://howlingpixel.com/wiki/Centripetal_Catmull%E2%80%93Rom_spline

https://stackoverflow.com/questions/9489736/catmull-rom-curve-with-no-cusps-and-no-self-intersections/23980479#23980479


----

[travis-img]: https://travis-ci.org/JeffreySarnoff/CentripetalCatmullRom.jl.svg?branch=master
[travis-url]: https://travis-ci.org/JeffreySarnoff/CentripetalCatmullRom.jl



[pkg-0.6-img]: http://pkg.julialang.org/badges/CentripetalCatmullRom_0.6.svg
[pkg-0.6-url]: http://pkg.julialang.org/?pkg=CentripetalCatmullRom&ver=0.6
[pkg-0.7-img]: http://pkg.julialang.org/badges/CentripetalCatmullRom_0.7.svg
[pkg-0.7-url]: http://pkg.julialang.org/?pkg=CentripetalCatmullRom&ver=0.7
