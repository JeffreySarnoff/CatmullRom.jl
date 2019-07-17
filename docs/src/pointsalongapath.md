
## Points along a path

A sequence of 2D, 3D .. nD points is required.  There is no limit on the number of coordinate dimensions.  
The first coordinate of each point become the abcissae (e.g. the `x` coordinate values).  The second \[, third etc.\]
become \[successive\] ordinates (e.g. the `ys`, `zs` ...).

Every point in a givne sequence must has the same number of constiuent coordinates.  Coordinates are considered to be values
along _orthonormal_ axes.  All ordinate axes are fitted with respect to the same abcissae. So, the arcs
that connect successive `y`s are arcs hewn from a succession of `(x_i, y_i)` ordered pairs and the arcs connecting successive
`z`s are arcs hewn from a succession of `(x_i, z_i)` ordered pairs.  It is easy to work with other axial pairings. To generate
arcs using the sequence of `(y_i, z_i)` pairs: `ys_zs = catmullrom( collect(zip(ys, zs)) )`.

----

The point sequence itself may be provided as a vector of points or as a tuple of points. 


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
