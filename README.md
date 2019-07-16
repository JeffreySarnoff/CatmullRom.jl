# CatmullRom.jl

### Centripetal parameterization for Catmull-Rom interpoint connectivity. 


#### Copyright ©&thinsp;2018-2019 by Jeffrey Sarnoff. &nbsp;&nbsp;  This work is released under The MIT License.


-----


[![Build Status](https://travis-ci.org/JeffreySarnoff/CatmullRom.jl.svg?branch=master)](https://travis-ci.org/JeffreySarnoff/CatmullRom.jl)
-----


  
## General Perspective

Catmull-Rom splines are a workhorse of computer graphics. Using the centripetal parameterization, they become a very handy general purpose tool for fast, attractive curvilinear blending. Often, they give interpoint "motion" a naturalistic feel.

----

<p align="center">
  <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/hair.PNG" width="750"> 
</p>

  <p align="center"><a href="http://www.cemyuksel.com/research/catmullrom_param">Cem Yuksel's Research</a></p>

----

<p align="center">
  <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/CR_Centripetal.png" width="400"> 
</p>

  <p align="center"><a href="http://www.cemyuksel.com/research/catmullrom_param/catmullrom_cad.pdf">Parameterization and Applications of Catmull-Rom Curves</a></p>


----

## Three functions are exported

### Open and Closed Curves

CatmullRom processes the extremal points of closed curves differently from open curves.
A curve in which the first and last points are identical is recognized as closed.
A function, `close_seq`, is available to ensure curves intended to be closed are made closed
in an exact and proper way. It is _good practice_ to use this function with closed curves,
and so assure they are crisp where the escribed path rejoins itself.
```
close_seq( points )            # this is the only function that may change some part of your data
                               # any change is limited to copying the first point into the last 
points = close_seq( points )   # (the same thing)
```

### Uniform Intermediation

There are two ways to connect path-adjacent points using Centripetal Catmull-Rom machinery.
The most often used iterplaces a given number of curvilinear waypoints between each adjacent
pair of original points.  All neighbors become connected by that given number of intermediating
places. Though the places differ, the proportional advancing between abcissae is consistent.
There is a default for this intermediation count, nonetheless trying a few different values
may help you visualize the significance that is of import.
```
crpoints = catmullrom( points )

crpoints = catmullrom( points, n_between_points )
```

### Arclength Relative Allocation

When the points' coordinates are spread differently along distinct axes, the interpoint
distances along one coordinate have a very different nature from the intercoordinate
spreads along another coordinate axis.  The distances separating adjacent point pairs
may vary substantively.  This is particularly true when working in higher dimensional
regions of an orthonormal coordinate space.  One may use more intermediating placements
between adjacent points that are relatively far apart, and fewer between adjacent points
that are in close relative proximity.

```
crpoints = catmullrom_by_arclength( points )

crpoints = catmullrom_by_arclength( points, (min_between_points, max_between_points) )
```

----

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




----

### points along a path

A sequence of 2D, 3D .. nD points is required.  There is no limit on the number of coordinate dimensions.  The sequence of values
given as the first coordinate of each point becomes the abcissae (the `x` coordinate values).  The second values become the
ordinates.  When there are more than two coordinates comprising each point, the second coordinate is as the `y` coordinate value
(or whatever coordinate is identified with the axis that follows e.g. by the right-hand rule).

Every point in a sequence of points has the same number of constiuent coordinates.  The first coordinate from each point,
in sequence, are the abcissae.  The remaining coordinates (e.g. the `y` or the `y` and `z` coordinates) are considered to be
values along orthonormal _ordinate_ axes.  All ordinate sequences are fitted with respect to the same abcissae. So, the arcs
that connect successive `y`s are arcs hewn from a succession of `(x_i, y_i)` ordered pairs and the arcs connecting successive
`z`s are arcs hewn from a succession of `(x_i, z_i)` ordered pairs.  It is easy to work with other axial pairings. To generate
arcs using the sequence of `(y_i, z_i)` pairs: `ys_zs = catmullrom( collect(zip(ys, zs)) )`.

The point sequence may be provided as a vector of points or as a tuple of points.  The points themselves may be vectors
of coordinate values or tuples of coordinate values.  While the points and their coordinates are manipulated internally,
that occurs without altering any values or sequences you use. We presuppose that you will have some carrier for each
`Point`, and a container for the points holding them as they are sequenced along whatever you deem their natural path.


----

|  Type used for a Point | example             |  coordinates are retrievable |  you support   | 
|:-----------------------|:--------------------|------------------------------|----------------|
|                        |                     |                                               |
|  small vector          | \[1.0, 3.5 \]       |   coord(point, i) = point[i] |   _builtin_    |
|  small tuple           | (1.0, 3.5)          |   coord(point, i) = point[i] |   _builtin_    |
|                        |                     |                                               |
|  StaticVector          | SVector( 1.0, 3.5 ) |   coord(point, i) = point[i] |   _builtin_    |
|  NamedTuple            | (x = 1.0, y = 3.5 ) |   coord(point, i) = point[i] |   _builtin_    |
|                        |                     |                                               |
|  struct                | Point(1.0, 3.5)     |   coord(point, i) = point[i] |   getindex     |
|                        |                     |                                               |

```
struct Point{T}
    x::T
    y::T
    z::T
end

function Base.getindex(point::Point{T}, i::Integer) where T
    if i == 1
       point.x
    elseif i == 2
       point.y
    elseif i == 3
       point.z
    else
       throw(DomainError("i must be 1, 2, or 3 (not $i)"))
    end
end
```



----

## Centripetal Catmull-Rom Examples <sup>[𝓪](#source)</sup>




  |                     shape                               |              detail                           |
  |:--------------------------------------------------------:|:-----------------------------------------------:|
  |                                                          |                                                 |
  | <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/CatmullRom_circle_dpihalf.png" width="300">  |      <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/CatmullRom_sectionofcircle.png" width="300">|
  |                                                          |                                                 |
   |    <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/teardrop_byarc.png" width="400">          |   <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/teardrop_section.PNG" width="250">         |
   |    <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/sphericalspiral.png" width="400">          |   <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/sphericalspiral_locus.png" width="130">    |  
  
  <p align="center"><sup><a name="source">𝓪</sup> <a href="https://github.com/JeffreySarnoff/CatmullRom.jl/tree/master/examples">using this package to generate some of these examples</a></p>  



----
### the first and last points are special

|    |   |
|:---------------------------------------------------------------------------------------------------------------------------|:--|
| <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/Catmull-Rom_Spline.png" width="500"> | [from Wikipedia](https://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline)  | 
| Catmull-Rom splines over two points are made with their neighbors. A new point preceeds your first and another follows your last. |
By appending points outside those given, the generated curve is bounded at your extremal points. |
| | This just happens with the internal flow, so if you add bounding points, use `catmullrom(points, n_between, extend=false)`.|

----

## hints
    
- If your points are disaggregated (e.g. all the `xs` in vec_of_xs, all the `ys` in vec_of_ys)
    - aggregate them this way `points = collect(zip(xs, ys, zs))`

- Often, abcissae (`xs`) are given in an ascending or in a descending sequence
    - `x[i-1] < x[i] < x[i+1]` or `x[i-1] > x[i] > x[i+1]`     

- With closed curves, expect one of these adjacency triplets   
    - `x[i-1] < x[i] > x[i+1]` or `x[i-1] > x[i] < x[i+1]`

----


### references

[Parameterization and Applications of Catmull-Rom Curves](http://www.cemyuksel.com/research/catmullrom_param/catmullrom_cad.pdf)

[The Centripetal Catmull-Rom Spline](https://howlingpixel.com/wiki/Centripetal_Catmull%E2%80%93Rom_spline)

[Catmull-Rom spline without cusps or self-intersections](https://stackoverflow.com/questions/9489736/catmull-rom-curve-with-no-cusps-and-no-self-intersections/23980479#23980479)

-----

[travis-img]: https://travis-ci.org/JeffreySarnoff/CatmullRom.jl.svg?branch=master
[travis-url]: https://travis-ci.org/JeffreySarnoff/CatmullRom.jl

[pkg-1.0-img]: http://pkg.julialang.org/badges/CatmullRom_1.0.svg
[pkg-1.0-url]: http://pkg.julialang.org/?pkg=CatmullRom&ver=1.0
