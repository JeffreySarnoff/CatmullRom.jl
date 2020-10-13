

parameter | represents        | value type  
:--------:|-------------------|-------------
   `ND`   | dimensions        |  `Signed`   
   `CT`   | coord type        |  `Number`   
   `KS`   | kind of spline    |  Singleton below   
   `PD`   | polynomial degree |  `Int`      
   
-----

_Related Symbols_   
```
const SplineKinds   = (:Bezier, :CatmullRom)
const SplineDegrees = (:Linear = 1, :Quadratic = 2, :Cubic = 3, :Quartic = 4)
```

