#  CentripetalCatmullRom.jl

### Centripetal parameterization for Catmull-Rom cubic interpolants. 


----

#### Copyright ©&thinsp;2018 by Jeffrey Sarnoff. &nbsp;&nbsp;  This work is released under The MIT License.


-----


[![][pkg-0.7-img]][pkg-0.7-url]  [![][travis-img]][travis-url]


-----

### Exports

## constructor

- catmullrom(points, interpolants)
    - `points` is a tuple of points-as-tuples
    - `interpolants` is a tuple of values from 0.0 to 1.0 (inclusive)
    -  yields interpolating points from points[2] through points[end-1] (inclusive)

## interpolation points

- `uniformsep(n)`
    - [0.0, b, c.., n-2, n-1, 1.0]
    - `c - b ~= (n-1) - (n-2)` 
    - _uniformly separated interpolants_

- `chebyshevsep(n)`
    - [0.0, b .. n-1, 1.0]
    - Chebyshev polynomial of the second kind, roots
    - roots of U(n), mapped into 0.0:1.0
    - these work better than roots T(n) with centripetal Catmull-Rom curves

## utilities

- `into01((xs...,))`, `into01([xs...,])`
    - maps values into 0.0:1.0, linearly
    - minimum(xs) --> 0.0, maximum(xs) --> 1.0

- `clamp01((xs...,))`, `clamp01([xs...,])`
    - forces values into 0.0..1.0

```julia
julia> using CentripetalCatmullRom

julia> interpolants = collect
julia> points2D = ([(sinpi(x),cospi(x)) for x=0.0f0:(0.25f0/3.0f0):0.25f0]...,)
((0.0f0, 1.0f0), (0.25881904f0, 0.9659258f0), (0.5f0, 0.8660254f0), (0.70710677f0, 0.70710677f0))

julia> polys = catmullrom_polys(points2D)
2-element Array{Poly{Float32},1}:
 Poly(0.25881904f0 + 0.25f0*x - 0.00060100853f0*x^2 - 0.008218035f0*x^3)   
 Poly(0.9659258f0 - 0.06698731f0*x - 0.03631732f0*x^2 + 0.0034040362f0*x^3)

julia> catmullrom(points2D,
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
