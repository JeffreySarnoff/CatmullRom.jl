

|  Centripetal Catmull-Rom Pathmaking   |
|:--------------------------------------|
|                                       |
|  `catmullrom( points )`               |
|                                       |
|  `catmullrom_by_arclength( points )`  |
|                                       |

## Uniform Intermediation

- catmullrom( points_along_a_path )
- catmullrom( points_along_a_path, n_arcs_between_neighbors )

There are two ways to connect path-adjacent points using Centripetal Catmull-Rom machinery.
The most often used iterplaces a given number of curvilinear waypoints between each adjacent
pair of original points.  All neighbors become connected by that given number of intermediating
places. Though the places differ, the proportional advancing between abcissae is consistent.
There is a default for this intermediation count, nonetheless trying a few different values
may help you visualize the significance that is of import.
```
crpoints = catmullrom( points )

crpoints = catmullrom( points, n_arcs_between_neighbors )
```
----

## Arclength Relative Allocation

When the points' coordinates are spread differently along distinct axes, the interpoint
distances along one coordinate have a very different nature from the intercoordinate
spreads along another coordinate axis.  The distances separating adjacent point pairs
may vary substantively.  This is particularly true when working in higher dimensional
regions of an orthonormal coordinate space.  One may use more intermediating placements
between adjacent points that are relatively far apart, and fewer between adjacent points
that are in close relative proximity.

```
crpoints = catmullrom_byarc( points )

crpoints = catmullrom_byarc( points, (min_arcs_between_points, max_arcs_between_points) )
```

----

```
using CatmullRom, Plots

result = catmullrom(points, n_arcs_per_pair)    # your points, how many new points to place between adjacents
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
 
result = catmullrom_byarc(points, (atleast_min_arcs_total, atmost_max_arcs_total))
                                  # min, max pertain to the whole path of the curve
xs, ys = result
plot(xs, ys)                      # or plot(result...)
```

----

## Open and Closed Curves

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

----
