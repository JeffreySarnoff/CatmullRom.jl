#=

| param  | represents | types used  |
|--------|------------|-------------|
| `ND`   | number of  |  `Signed`   |
|        | dimensions |             |
|--------|------------|-------------|
| `CT`   | coordinate |  `Number`   |
|        | type       |             |
|--------|------------|-------------|
| `KS`   | kind of    | :Bezier     |
|        | spline     | :CatmullRom |
|--------|------------|-------------|
| `PD`   | polynomial |  `Int`      |
|        | degree     |             |
|--------|------------|-------------|


=#
struct SplineSpan{ND, CT
struct QuadraticBezier{S,N,T} <: AbstractCurve{S,N,T} end
struct
function approx_arclength(::Type{Bezier), ::Type{Quadratic},
#=
def approximateQuadraticArcLengthC(pt1, pt2, pt3):
    """Calculates the arc length for a quadratic Bezier segment.
    Uses Gauss-Legendre quadrature for a branch-free approximation.
    See :func:`calcQuadraticArcLength` for a slower but more accurate result.
    Args:
        pt1: Start point of the Bezier as a complex number.
        pt2: Handle point of the Bezier as a complex number.
        pt3: End point of the Bezier as a complex number.
    Returns:
        Approximate arc length value.
    """
    # This, essentially, approximates the length-of-derivative function
    # to be integrated with the best-matching fifth-degree polynomial
    # approximation of it.
    #
    #https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss.E2.80.93Legendre_quadrature
    # abs(BezierCurveC[2].diff(t).subs({t:T})) for T in sorted(.5, .5±sqrt(3/5)/2),
    # weighted 5/18, 8/18, 5/18 respectively.
    v0 = abs(-0.492943519233745*pt1 + 0.430331482911935*pt2 + 0.0626120363218102*pt3)
    v1 = abs(pt3-pt1)*0.4444444444444444
    v2 = abs(-0.0626120363218102*pt1 - 0.430331482911935*pt2 + 0.492943519233745*pt3)
    return v0 + v1 + v2
    
    =#
