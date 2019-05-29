"""
    CatmullRom

Interpolate using the centripetal parameterization of Catmull-Rom splines.
This is the form that does not crimp, provides taughtness, and looks good.
"""
module CatmullRom

export catmullrom,    # populates with points placed between those given
       extendbounds   # extrapolate just beyond endpoints from adjacent points

# available as `CatmullRom.extendbounds`, `CatmullRom.thiele3`, etc
# reflect, linear, quadratic, thiele3, thiele4  # also for interpolation

include("centripetal/catmullrom.jl")
include("centripetal/support.jl")
include("centripetal/arcbased.jl")

include("fewpoints/twopoints.jl")
include("fewpoints/threepoints.jl")
include("fewpoints/fourpoints.jl")

end # CatmullRom
