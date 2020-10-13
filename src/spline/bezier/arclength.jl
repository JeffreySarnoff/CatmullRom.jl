# from ob/master/Lib/fontTools/misc/bezierTools.py

function estimate_bezier_arclength(r0::T, r1::T, r2::T, r3::T) where {T<:Real}
    c0, c1, c2, c3 = complex.((r0, r1, r2, r3))
    return estimate_bezier_arclength(c0, c1, c2, c3)
end

function estimate_bezier_arclength(c0::T, c1::T, c2::T, c3::T) where {T<:Complex}
    v0 = abs(c1-c0) * 0.15
    v1 = abs(-0.558983582205757*c0 + 0.325650248872424*c1 + 0.208983582205757*c2 + 0.024349751127576*c3)
    v2 = abs(c3-c0+c2-c1) * 0.26666666666666666
    v3 = abs(-0.024349751127576*c0 - 0.208983582205757*c1 - 0.325650248872424*c2 + 0.558983582205757*c3)
    v4 = abs(c3-c2) * 0.15
    return real(sum((v0, v1, v2, v3, v4)))
end    

function estimate_bezier_arclength(r0::T, r1::T, r2::T, r3::T) where {T<:Real}
    c0, c1, c2, c3 = complex.((r0, r1, r2, r3))
    return estimate_bezier_arclength(c0, c1, c2, c3)
end

function estimate_bezier_arclength(c0::T, c1::T, c2::T, c3::T) where {T<:Complex}
    v0 = abs(c1-c0) * 0.15
    v1 = abs(-0.558983582205757*c0 + 0.325650248872424*c1 + 0.208983582205757*c2 + 0.024349751127576*c3)
    v2 = abs(c3-c0+c2-c1) * 0.26666666666666666
    v3 = abs(-0.024349751127576*c0 - 0.208983582205757*c1 - 0.325650248872424*c2 + 0.558983582205757*c3)
    v4 = abs(c3-c2) * 0.15
    return real(sum((v0, v1, v2, v3, v4)))
end    


#=



def calcCubicArcLength(pt1, pt2, pt3, pt4, tolerance=0.005):
    """Calculates the arc length for a cubic Bezier segment.
    Whereas :func:`approximateCubicArcLength` approximates the length, this
    function calculates it by "measuring", recursively dividing the curve
    until the divided segments are shorter than ``tolerance``.
    Args:
        pt1,pt2,pt3,pt4: Control points of the Bezier as 2D tuples.
        tolerance: Controls the precision of the calcuation.
    Returns:
        Arc length value.
    """
    return calcCubicArcLengthC(complex(*pt1), complex(*pt2), complex(*pt3), complex(*pt4), tolerance)


def _split_cubic_into_two(p0, p1, p2, p3):
    mid = (p0 + 3 * (p1 + p2) + p3) * .125
    deriv3 = (p3 + p2 - p1 - p0) * .125
    return ((p0, (p0 + p1) * .5, mid - deriv3, mid),
            (mid, mid + deriv3, (p2 + p3) * .5, p3))

def _calcCubicArcLengthCRecurse(mult, p0, p1, p2, p3):
	arch = abs(p0-p3)
	box = abs(p0-p1) + abs(p1-p2) + abs(p2-p3)
	if arch * mult >= box:
		return (arch + box) * .5
	else:
		one,two = _split_cubic_into_two(p0,p1,p2,p3)
		return _calcCubicArcLengthCRecurse(mult, *one) + _calcCubicArcLengthCRecurse(mult, *two)

def calcCubicArcLengthC(pt1, pt2, pt3, pt4, tolerance=0.005):
    """Calculates the arc length for a cubic Bezier segment.
    Args:
        pt1,pt2,pt3,pt4: Control points of the Bezier as complex numbers.
        tolerance: Controls the precision of the calcuation.
    Returns:
        Arc length value.
    """
    mult = 1. + 1.5 * tolerance # The 1.5 is a empirical hack; no math
    return _calcCubicArcLengthCRecurse(mult, pt1, pt2, pt3, pt4)


epsilonDigits = 6
epsilon = 1e-10


def _dot(v1, v2):
    return (v1 * v2.conjugate()).real


def _intSecAtan(x):
    # In : sympy.integrate(sp.sec(sp.atan(x)))
    # Out: x*sqrt(x**2 + 1)/2 + asinh(x)/2
    return x * math.sqrt(x**2 + 1)/2 + math.asinh(x)/2


def calcQuadraticArcLength(pt1, pt2, pt3):
    """Calculates the arc length for a quadratic Bezier segment.
    Args:
        pt1: Start point of the Bezier as 2D tuple.
        pt2: Handle point of the Bezier as 2D tuple.
        pt3: End point of the Bezier as 2D tuple.
    Returns:
        Arc length value.
    Example::
        >>> calcQuadraticArcLength((0, 0), (0, 0), (0, 0)) # empty segment
        0.0
        >>> calcQuadraticArcLength((0, 0), (50, 0), (80, 0)) # collinear points
        80.0
        >>> calcQuadraticArcLength((0, 0), (0, 50), (0, 80)) # collinear points vertical
        80.0
        >>> calcQuadraticArcLength((0, 0), (50, 20), (100, 40)) # collinear points
        107.70329614269008
        >>> calcQuadraticArcLength((0, 0), (0, 100), (100, 0))
        154.02976155645263
        >>> calcQuadraticArcLength((0, 0), (0, 50), (100, 0))
        120.21581243984076
        >>> calcQuadraticArcLength((0, 0), (50, -10), (80, 50))
        102.53273816445825
        >>> calcQuadraticArcLength((0, 0), (40, 0), (-40, 0)) # collinear points, control point outside
        66.66666666666667
        >>> calcQuadraticArcLength((0, 0), (40, 0), (0, 0)) # collinear points, looping back
        40.0
    """
    return calcQuadraticArcLengthC(complex(*pt1), complex(*pt2), complex(*pt3))


def calcQuadraticArcLengthC(pt1, pt2, pt3):
    """Calculates the arc length for a quadratic Bezier segment.
    Args:
        pt1: Start point of the Bezier as a complex number.
        pt2: Handle point of the Bezier as a complex number.
        pt3: End point of the Bezier as a complex number.
    Returns:
        Arc length value.
    """
    # Analytical solution to the length of a quadratic bezier.
    # I'll explain how I arrived at this later.
    d0 = pt2 - pt1
    d1 = pt3 - pt2
    d = d1 - d0
    n = d * 1j
    scale = abs(n)
    if scale == 0.:
        return abs(pt3-pt1)
    origDist = _dot(n,d0)
    if abs(origDist) < epsilon:
        if _dot(d0,d1) >= 0:
            return abs(pt3-pt1)
        a, b = abs(d0), abs(d1)
        return (a*a + b*b) / (a+b)
    x0 = _dot(d,d0) / origDist
    x1 = _dot(d,d1) / origDist
    Len = abs(2 * (_intSecAtan(x1) - _intSecAtan(x0)) * origDist / (scale * (x1 - x0)))
    return Len


def approximateQuadraticArcLength(pt1, pt2, pt3):
    """Calculates the arc length for a quadratic Bezier segment.
    Uses Gauss-Legendre quadrature for a branch-free approximation.
    See :func:`calcQuadraticArcLength` for a slower but more accurate result.
    Args:
        pt1: Start point of the Bezier as 2D tuple.
        pt2: Handle point of the Bezier as 2D tuple.
        pt3: End point of the Bezier as 2D tuple.
    Returns:
        Approximate arc length value.
    """
    return approximateQuadraticArcLengthC(complex(*pt1), complex(*pt2), complex(*pt3))


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


def calcQuadraticBounds(pt1, pt2, pt3):
    """Calculates the bounding rectangle for a quadratic Bezier segment.
    Args:
        pt1: Start point of the Bezier as a 2D tuple.
        pt2: Handle point of the Bezier as a 2D tuple.
        pt3: End point of the Bezier as a 2D tuple.
    Returns:
        A four-item tuple representing the bounding rectangle ``(xMin, yMin, xMax, yMax)``.
    Example::
        >>> calcQuadraticBounds((0, 0), (50, 100), (100, 0))
        (0, 0, 100, 50.0)
        >>> calcQuadraticBounds((0, 0), (100, 0), (100, 100))
        (0.0, 0.0, 100, 100)
    """
    (ax, ay), (bx, by), (cx, cy) = calcQuadraticParameters(pt1, pt2, pt3)
    ax2 = ax*2.0
    ay2 = ay*2.0
    roots = []
    if ax2 != 0:
        roots.append(-bx/ax2)
    if ay2 != 0:
        roots.append(-by/ay2)
    points = [(ax*t*t + bx*t + cx, ay*t*t + by*t + cy) for t in roots if 0 <= t < 1] + [pt1, pt3]
    return calcBounds(points)


def approximateCubicArcLength(pt1, pt2, pt3, pt4):
    """Approximates the arc length for a cubic Bezier segment.
    Uses Gauss-Lobatto quadrature with n=5 points to approximate arc length.
    See :func:`calcCubicArcLength` for a slower but more accurate result.
    Args:
        pt1,pt2,pt3,pt4: Control points of the Bezier as 2D tuples.
    Returns:
        Arc length value.
    Example::
        >>> approximateCubicArcLength((0, 0), (25, 100), (75, 100), (100, 0))
        190.04332968932817
        >>> approximateCubicArcLength((0, 0), (50, 0), (100, 50), (100, 100))
        154.8852074945903
        >>> approximateCubicArcLength((0, 0), (50, 0), (100, 0), (150, 0)) # line; exact result should be 150.
        149.99999999999991
        >>> approximateCubicArcLength((0, 0), (50, 0), (100, 0), (-50, 0)) # cusp; exact result should be 150.
        136.9267662156362
        >>> approximateCubicArcLength((0, 0), (50, 0), (100, -50), (-50, 0)) # cusp
        154.80848416537057
    """
    return approximateCubicArcLengthC(complex(*pt1), complex(*pt2), complex(*pt3), complex(*pt4))


def approximateCubicArcLengthC(pt1, pt2, pt3, pt4):
    """Approximates the arc length for a cubic Bezier segment.
    Args:
        pt1,pt2,pt3,pt4: Control points of the Bezier as complex numbers.
    Returns:
        Arc length value.
    """
    # This, essentially, approximates the length-of-derivative function
    # to be integrated with the best-matching seventh-degree polynomial
    # approximation of it.
    #
    # https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss.E2.80.93Lobatto_rules

    # abs(BezierCurveC[3].diff(t).subs({t:T})) for T in sorted(0, .5±(3/7)**.5/2, .5, 1),
    # weighted 1/20, 49/180, 32/90, 49/180, 1/20 respectively.
    v0 = abs(pt2-pt1)*.15
    v1 = abs(-0.558983582205757*pt1 + 0.325650248872424*pt2 + 0.208983582205757*pt3 + 0.024349751127576*pt4)
    v2 = abs(pt4-pt1+pt3-pt2)*0.26666666666666666
    v3 = abs(-0.024349751127576*pt1 - 0.208983582205757*pt2 - 0.325650248872424*pt3 + 0.558983582205757*pt4)
    v4 = abs(pt4-pt3)*.15

    return v0 + v1 + v2 + v3 + v4


def calcCubicBounds(pt1, pt2, pt3, pt4):
    """Calculates the bounding rectangle for a quadratic Bezier segment.
    Args:
        pt1,pt2,pt3,pt4: Control points of the Bezier as 2D tuples.
    Returns:
        A four-item tuple representing the bounding rectangle ``(xMin, yMin, xMax, yMax)``.
    Example::
        >>> calcCubicBounds((0, 0), (25, 100), (75, 100), (100, 0))
        (0, 0, 100, 75.0)
        >>> calcCubicBounds((0, 0), (50, 0), (100, 50), (100, 100))
        (0.0, 0.0, 100, 100)
        >>> print("%f %f %f %f" % calcCubicBounds((50, 0), (0, 100), (100, 100), (50, 0)))
        35.566243 0.000000 64.433757 75.000000
    """
    (ax, ay), (bx, by), (cx, cy), (dx, dy) = calcCubicParameters(pt1, pt2, pt3, pt4)
    # calc first derivative
    ax3 = ax * 3.0
    ay3 = ay * 3.0
    bx2 = bx * 2.0
    by2 = by * 2.0
    xRoots = [t for t in solveQuadratic(ax3, bx2, cx) if 0 <= t < 1]
    yRoots = [t for t in solveQuadratic(ay3, by2, cy) if 0 <= t < 1]
    roots = xRoots + yRoots

    points = [(ax*t*t*t + bx*t*t + cx * t + dx, ay*t*t*t + by*t*t + cy * t + dy) for t in roots] + [pt1, pt4]
    return calcBounds(points)

=#
