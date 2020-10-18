"""
    fast_catmullrom_arclength(pt0, pt1, pt2, pt3)

Given 4 ND points along a centripetal Catmull-Rom span,
roughly approximate the arclength of the curvilinear segment
determined by the two bounding points and their tangents
[the arc between pt1 and pt2].
This well-behaved approximation was developed by Jens Gravesen
```
    (2*corddist + (n-1)*bezdist)/(n+1), n is degree of the curve
    deg=2 --> (2*corddist + bezdist)/(3)
    deg=3 --> (2*corddist + 2*bezdist)/(4) --> (corddist + bezdist)/2
    deg=4 --> (2*corddidt + 3*bezdist)/(5)
```
"""
function fast_catmullrom_arclength(p0::T, p1::T, p2::T, p3::T) where T
    b0, b1, b2, b3 = catmullrom_as_bezier(p0, p1, p2, p3)
    cordal_dist = norm(b3 .- b0)
    bezier_dist = norm(b1 .- b0) + norm(b2 .- b1) + norm(b3 .- b2)
    
    return (cordal_dist + bezier_dist) * 0.5
end

"""
    approx_catmullrom_arclength(pt0, pt1, pt2, pt3)

Given 4 ND points along a centripetal Catmull-Rom span,
approximate the arclength of the curvilinear segment
determined by the two bounding points and their tangents
[the arc between pt1 and pt2].

This, essentially, approximates the length-of-derivative function
to be integrated with the best-matching seventh-degree polynomial
approximation of it.
https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss.E2.80.93Lobatto_rules
   for t in [0, sqrt(0.5-3/7)/2, sqrt(0.5+3/7)/2, .5, 1],
   weighted 1/20, 49/180, 32/90, 49/180, 1/20 respectively.
"""
function approx_catmullrom_arclength(p1::T, p2::T, p3::T, p4::T) where T
    b1, b2, b3, b4 = catmullrom_as_bezier(p1, p2, p3, p4)
    c1 = complex(b1...)
    c2 = complex(b2...)
    c3 = complex(b3...)
    c4 = complex(b4...)
    v0 = abs(c2-c1) * 0.15
    v1 = abs(-0.558983582205757*c1 +
              0.325650248872424*c2 +
              0.208983582205757*c3 +
              0.024349751127576*c4)
    v2 = abs(c4-c1+c3-c2) * 0.26666666666666666
    v3 = abs(-0.024349751127576*c1 -
              0.208983582205757*c2 -
              0.325650248872424*c3 +
              0.558983582205757*c4)
    v4 = abs(c4-c3) * 0.15

    return v0 + v1 + v2 + v3 + v4
end

function catmullrom_as_bezier(p0::T, p1::T, p2::T, p3::T) where T
    b0 = p1
    b3 = p2
    
    d1 = norm(p1 .- p0); d1a = sqrt(d1)
    d2 = norm(p2 .- p1); d2a = sqrt(d2)
    d3 = norm(p3 .- p2); d3a = sqrt(d3)
    
    b1n = @. (d1 * p2) - (d2 * p0) + ((2*d1 + 3*d1a*d2a+d2) * p1)
    b1d = 3*d1a*(d1a+d2a)
    b1 = b1n ./ b1d
    b2n = @. (d3 * p1) - (d2 * p3) + ((2*d3 + 3*d3a*d2a+d2) * p2)
    b2d = 3*d3a*(d3a+d2a)
    b2 = b2n ./ b2d
    
    return b0, b1, b2, b3
end
