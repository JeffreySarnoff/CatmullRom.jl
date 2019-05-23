"""
    CatmullRom

This package provides open and closed centripetal Catmull-Rom splines.
The centripetal parameterization of Catmull-Rom splines works best in use:
it does not crimp and provides a useful amount of relative taughtness.
"""
module CatmullRom

export catmullrom,
       # extrapolate just beyond endpoints from 2 or 3 given, adjacent points
       reflect, linear, quadratic, thiele3, thiele4  # also for interpolation
       
include("fewpoints/twopoints.jl")
include("fewpoints/threepoints.jl")
include("fewpoints/fourpoints.jl")


end # CatmullRom
