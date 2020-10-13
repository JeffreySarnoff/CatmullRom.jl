

parameter | represents        | value type  
:--------:|-------------------|-------------
   `ND`   | dimensions        |  `Signed`   
   `CT`   | coord type        |  `Number`   
   `KS`   | kind of spline    | spline symbols  
                                :CatmullRom)     
   `PD`   | polynomial degree |  `Int`      
   
   
const SplineKinds = (:Bezier, :CatmullRom)
const SplineFits  = (:Linear, :Quadratic, :Cubic, :Quartic)

