# CatmullRom.jl

### Centripetal parameterization for Catmull-Rom interpoint connectivity. 


#### Copyright ©&thinsp;2018-2019 by Jeffrey Sarnoff. &nbsp;&nbsp;  This work is released under The MIT License.


-----


[![Build Status](https://travis-ci.org/JeffreySarnoff/CatmullRom.jl.svg?branch=master)](https://travis-ci.org/JeffreySarnoff/CatmullRom.jl)
-----


  
## General Perspective

Catmull-Rom splines are a workhorse of computer graphics and a very handy general purpose tool for fast, attractive blending.  Using the centripetal parameterization gives a more naturalistic feel to the interpoint "motion".

|  The Centripetal Catmull-Rom curve is Blue |
|:-----------------------------------------------------------------------------------------------------------------------:|
|  <img src="https://github.com/JeffreySarnoff/CatmullRom.jl/blob/master/examples/assets/CR_Centripetal.png" width="300"> |
|                                                                                                                         |
| from "Parameterization and Applications of Catmull-Rom Curves" Cem Yuksel, Scott Schaefer, John Keyser  | 

## Use



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
