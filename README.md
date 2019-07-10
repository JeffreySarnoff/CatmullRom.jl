# CatmullRom.jl

### Centripetal parameterization for Catmull-Rom interpoint connectivity. 


#### Copyright ©&thinsp;2018-2019 by Jeffrey Sarnoff. &nbsp;&nbsp;  This work is released under The MIT License.


-----


[![Build Status](https://travis-ci.org/JeffreySarnoff/CatmullRom.jl.svg?branch=master)](https://travis-ci.org/JeffreySarnoff/CatmullRom.jl)
-----


  
## General Perspective

Catmull-Rom splines are a workhorse of computer graphics. Using the centripetal parameterization, they become a very handy general purpose tool for fast, attractive curvilinear blending. Often, they give interpoint "motion" a naturalistic feel.

|  The Centripetal Catmull-Rom Curve is Blue |
|:-----------------------------------------------------------------------------------------------------------------------:|
|  <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/CR_Centripetal.png" width="300"> |
|                                                                                                                         |
| [Parameterization and Applications of Catmull-Rom Curves by Cem Yuksel, Scott Schaefer, John Keyser](http://www.cemyuksel.com/research/catmullrom_param/catmullrom_cad.pdf)  | 

## Use

A sequence of 2D, 3D .. nD points is required.  There is no limit on the number of coordinate dimensions. The point sequence may be provided as a vector of points or as a tuple of points.  The points themselves may be vectors of coordinate values or tuples of coordinate values.


| <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/Catmull-Rom_Spline.png" width="300">  |
|:---------------------------------------------------------------------------------------------------------------------------:|
| [image from Wikipedia entry on Catmull-Rom Splines](https://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline) |

The point sequence that is used will be modified in place (a point will be prepended and a point will be postpended).  This is necessary because each Catmull-Rom spline uses four points to characterize the arc connecting the center two of those points. By appending a new first and new last point, the initial and final points of your sequence become the initial and final points in the result.


```
using CatmullRom, Plots

result = catmullrom(points, n_inbetween_points)  # your points, how many new points to place between adjacents
                                                 # result is a vector of coordinates, e.g. [xs, ys, zs]
plot(result...,)
```

When your points have nonuniform separation, or separation extents vary with coordinate dimension,
it is of benefit to allocate more of the new inbetween points where there are relatively greater
distances between your adjacent points.  The most appropriate measure for comparison and weighting
is interpoint arclength.  This package implements a well-behaved approximation to Catmull-Rom
arclengths appropriate to the centripetal parameterization.  You can use this directly.

```
using CatmullRom, Plots

result = catmullrom_byarc(points) # result is a vector of coordinates, e.g. [xs, ys, zs]
 
result = catmullrom_byarc(points, arcpoints_min=at_least, arcpoints_max=at_most)
                                  # specify the range of inbetween points used
xs, ys = result
plot(xs, ys)                      # plot(result...,)
```

If your points exist as separate coordinate vectors, aggregate them this way
`points = collect(zip(xs, ys, zs))`


## Centripetal Catmull-Rom Examples


  |                     shape                               |              detail                           |
  |:--------------------------------------------------------:|:-----------------------------------------------:|
  |                                                          |                                                 |
  | <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/CatmullRom_circle_dpihalf.png" width="300">  |      <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/CatmullRom_sectionofcircle.png" width="300">|
  |                                                          |                                                 |
   |    <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/teardrop_byarc.png" width="400">          |   <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/teardrop_section.PNG" width="250">         |
   |    <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/sphericalspiral.png" width="400">          |   <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/sphericalspiral_locus.png" width="130">    |  
  
  

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
