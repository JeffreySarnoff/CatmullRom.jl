

   parameter  | represents         | value type  
:------------:|--------------------|-------------
   __`ND`__   | dimensionalty      |  `Signed`   
   __`CT`__   | typeof(coordinate) |  `Number`   
   __`PD`__   | polynomial degree  |  `Int`      
   __`KS`__   | kind of spline     |  singelton   


```

----

#### these Kinds of Splines are internally supported


```
#=
    All Splines share the AbstractSpline abstraction.
    All Linear Splines share the AkoLinearSpline abstraction
    
    An unformed  spline does not join adjcent points with any algebric-geometric entity.
      An unformed spline does contain the sequence of coordinate-valued pointss
         which were given in its construction and knows whether it is cyclically
         connected (last point becomes the tail for a segment that reaches the first point)
         or that the initial and the final points are disconnected, separated.         
         
      Unformed splines are constructed from
          - a finite set of points within a shared space
              - replicated points are permissible, though unusual
          - an associatable total ordering over the unique points
          - the points gathered (or otherwise referencible) 
               in the sequence implicit in the ordering
               (  identical points are points with identical coordinate fields,
                  each point's corresponding fields holding eqivalent values)
                  identical points maitain thier relative order
                  as originally given to the constructor --
                  if their was no discernable relative order,
                  the first availble becomes the predecessor to the next. )
               
    A linear    spline uses line segments to join adjacent points.
    A quadratic spline uses second degree curve segments to join adjacent points.
    A cubic     spline uses segments of cubic polynomials to join adjacent points.
    A quartic   spline uses segments of 4th degree polynomials to join adjacent points.
    A quintic   spline uses segments of 5th degree polynomials to join adjacent points.
                

- the points made explicity available in that sequence.
           
           sequence of points that obtians
          
        Each point holds its constuent coordinate's fields. 
        Each constitutive coordianate is well-valued.
        
        Every constituent coordinate of each point is well-valued
        ach point's constiuent cooordinates are well-valued
        speifically valued
        values  with its coordinates -with-valued-coore onto which
        and an associated total ordering which unambiguously provides them sequenced.
        ssoociates is possibly circurlarelaborated    
      points that are adjacent in their given sequence
    
=#
```

```
```
abstract type AbstractSpline{ND,PD} end`

abstract type AkoLinearSpline{ND,PD} <: AbstractSpline{ND,PD} end
```
```

#=
    _Where you see "Ako", read it as "A kind of"._
    "AkoLinearSpline" is read "A kind of LinearSpline".
    
=#
abstract type AkoOrderedSpline{ND,PD}   <: AbstractSpline{ND,PD} end
abstract type AkoLinearSpline{ND,PD}    <: AbstractSpline{ND,PD} end
abstract type AkoQuadraticSpline{ND,PD} <: AbstractSpline{ND,PD} end
abstract type AkoCubicSpline{ND,PD}     <: AbstractSpline{ND,PD} end
abstract type AkoQuarticSpline{ND,PD}   <: AbstractSpline{ND,PD} end
abstract type AkoQuinticSpline{ND,PD}   <: AbstractSpline{ND,PD} end
```
-----

_Related Symbols_   

```julia
const SplineKinds   = (:Bezier, :CatmullRom)0
const SplineDegrees = (:Linear = 1, :Quadratic = 2, :Cubic = 3, :Quartic = 4)
```

