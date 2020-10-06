## What works as a `Point`?

#### Points reside in your analytic space

- How many dimensions are available in this space?
    - As many as there are distinct coordinate values contained in a `Point`.
    - a point that carries two values, e.g. `_x_, _y_`, is used within
      a 2D space or within a 2D projection of a larger dimensioned space
    - a point that carries three values is used within a 3D space
      (or a 3D projection from a larger dimensioned space)
      
- CatmullRom functions are designed with the assumption that each coordinate
  maps along a coordinate axis and these axes are algebraically distinct,
  geometrically orthogonal.
  
 - There is an initial axis, that to which the first coordinate of a point maps.
     - The first coordinate of all points (during a call) become abcissae.
     - The remaining coordinate[s] is[are] processed as the ordinate[s].
     
 - Working with domains that span 2 or more axes is not supported by this version.
 
## How are multidimensional ordinates processed?

### Interpolating a 1D sequence

avec = [ a₁ a₂ a₃ a₄ .. aᵢ₋₃ aᵢ₋₂ aᵢ₋₁ aᵢ aᵢ₊₁ aᵢ₊₂ aᵢ₊₃ .. aₙ₋₄ aₙ₋₃ aₙ₋₂ aₙ₋₁ aₙ ]

### Interpolating a 2D sequence

### Interpolating a 3D sequence

### Interpolating a nD sequence

 
