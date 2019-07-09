# CatmullRom.jl

### Centripetal parameterization for Catmull-Rom interpoint connectivity. 


#### Copyright ©&thinsp;2018-2019 by Jeffrey Sarnoff. &nbsp;&nbsp;  This work is released under The MIT License.


-----


[![Build Status](https://travis-ci.org/JeffreySarnoff/CatmullRom.jl.svg?branch=master)](https://travis-ci.org/JeffreySarnoff/CatmullRom.jl)
-----


## Centripetal Catmull-Rom Interpolation


  |                     shape                               |              detail                           |
  |:--------------------------------------------------------:|:-----------------------------------------------:|
  |                                                          |                                                 |
  | <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/CatmullRom_circle_dpihalf.png" width="300">  |      <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/CatmullRom_sectionofcircle.png" width="300">|
  |                                                          |                                                 |
   |    <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/teardrop_byarc.png" width="400">          |   <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/teardrop_section.PNG" width="250">         |      |    <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/sphericalspiral.png" width="400">          |   <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/sphericalspiral_locus.png" width="130">                                         |
  
  


### Use

```
using CatmullRom

xs = [..]  # these are the abcissae
ys = [..]  # there is at least one ordinate sequence
zs = [..}  # and at most, as many as you choose to use
           # each ordinate dimension is independent of the others
           # all ordinate sequences match to the abcissa sequence

# If your n-d curve is intended to be an open curve
# make sure that the first and last abcissae values are different
#
# If your n-d curve is intended to be a closed curve
# make sure that the first and last abcissae values are identical
# make sure that the first and last ordinate values for each ordinate dimension are identical

# each ordinate sequence must have the same number of values as there are abcissae
@assert length(xs) == length(ys) [ == length(zs) ]

n_interpoint_arcs = 42   # up to you, try different values to see what is best given the context
                         # use even numbers (so >= 2) and 12, 16, 24, 50, 64, 120 might work
                         
thepoints = zip(xs, ys [, zs ...])  # with (xs, ys) or (xs, ys, zs) you do not need to zip

crpoints = catmullrom(points, n_interpoint_arcs)                         
crpoints = catmullrom(xs, ys, n_interpoint_arcs)                         
crpoints = catmullrom(xs, ys, zs, n_interpoint_arcs)                         

# crpoints is a vector of vectors, on for each coordinate dimension
# crpoints includes the original points and adds all the interpoint placements

```


```julia

julia> points = ([(sinpi(x),cospi(x)) for x=0.0f0:(0.25f0/3.0f0):0.25f0]...,)
((0.0f0, 1.0f0), (0.25881904f0, 0.9659258f0), (0.5f0, 0.8660254f0), (0.70710677f0, 0.70710677f0))

julia> xys = CatmullRom.points_to_coords(points)
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
