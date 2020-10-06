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

```julia
xs = [x₁ x₂ x₃ x₄ x₅ .. xₖ₋₃ xₖ₋₂ xₖ₋₁ xₖ xₖ₊₁ xₖ₊₂ xₖ₊₃ ..  xₙ₋₄ xₙ₋₃ xₙ₋₂ xₙ₋₁ xₙ][:]
ys = [y₁ y₂ y₃ y₄ y₅ .. yₖ₋₃ yₖ₋₂ yₖ₋₁ yₖ yₖ₊₁ yₖ₊₂ yₖ₊₃ ..  yₙ₋₄ yₙ₋₃ yₙ₋₂ yₙ₋₁ yₙ][:]
zs = [z₁ z₂ z₃ z₄ z₅ .. zₖ₋₃ zₖ₋₂ zₖ₋₁ zₖ zₖ₊₁ zₖ₊₂ zₖ₊₃ ..  zₙ₋₄ zₙ₋₃ zₙ₋₂ zₙ₋₁ zₙ][:]
ws = [w₁ w₂ w₃ w₄ w₅ .. wₖ₋₃ wₖ₋₂ wₖ₋₁ wₖ wₖ₊₁ wₖ₊₂ wₖ₊₃ ..  wₙ₋₄ wₙ₋₃ wₙ₋₂ wₙ₋₁ wₙ][:]

dx = diff(xs) # xs[2:end] .- xs[1:end-1]
dy = diff(ys)
dz = diff(zs)
dw = diff(ws)

# chordal interpoint distances
dxy  = dist((x1,y1), (x2,y2)), dist((x2,y2), (x3,y3)) ...
dxz  = dist((x1,z1), (x2,z2)), dist((x2,z2), (x3,z3)) ...
dyz  = dist((y1,z1), (y2,z2)), dist((y2,z2), (y3,z3)) ...
dxyz = dist((x1,y1,z1), (x2,y2,z2)), dist((x2,y2,z2), (x3,y3,z3)) ...

myx = dy ./ dx
mzx = dz ./ dx
mwx = dw ./ dx
mzy = dz ./ dy
mwy = dw ./ dy
mwz = dw ./ dz

hs = [h₁ h₂ h₃ h₄ h₅ .. hₖ₋₃ hₖ₋₂ hₖ₋₁ hₖ hₖ₊₁ hₖ₊₂ hₖ₊₃ ..  hₙ₋₄ hₙ₋₃ hₙ₋₂ hₙ₋₁ hₙ][:]
ms = [m₁ m₂ m₃ m₄ m₅ .. mₖ₋₃ mₖ₋₂ mₖ₋₁ mₖ mₖ₊₁ mₖ₊₂ mₖ₊₃ ..  mₙ₋₄ mₙ₋₃ mₙ₋₂ mₙ₋₁ mₙ][:]
rs = [r₁ r₂ r₃ r₄ r₅ .. rₖ₋₃ rₖ₋₂ rₖ₋₁ rₖ rₖ₊₁ rₖ₊₂ rₖ₊₃ ..  rₙ₋₄ rₙ₋₃ rₙ₋₂ rₙ₋₁ rₙ][:]
ts = [t₁ t₂ t₃ t₄ t₅ .. tₖ₋₃ tₖ₋₂ tₖ₋₁ tₖ tₖ₊₁ tₖ₊₂ tₖ₊₃ ..  tₙ₋₄ tₙ₋₃ tₙ₋₂ tₙ₋₁ tₙ][:]
```

### Interpolating a 1D sequence

### Interpolating a 2D sequence

### Interpolating a 3D sequence

### Interpolating a nD sequence

 
