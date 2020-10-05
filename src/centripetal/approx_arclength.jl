"""
    approx_catmullrom_arclength(pt0, pt1, pt2, pt3)

Given 4 sequential points, approximate the arclength between pt1 and pt2.

This approximates the length of curvilinear segment determined by
the two bounding points [pt1, pt2] and the tangents determined by
the two bounding point pairs [(pt0, pt1), (pt2, pt3)].

This well-behaved approximation was developed by Jens Gravesen
```
    (2*corddist + (n-1)*bezdist)/(n+1), n is degree of the curve
    deg=2 --> (2*corddist + bezdist)/(3)
    deg=3 --> (2*corddist + 2*bezdist)/(4) --> (corddist + bezdist)/2
    deg=4 --> (2*corddsit + 3*bezdist)/(5)
```
"""
function approx_catmullrom_arclength(p0::T, p1::T, p2::T, p3::T) where T
    b0, b1, b2, b3 = catmullrom_as_bezier(p0, p1, p2, p3)
    cordal_dist = norm(b3 .- b0)
    bezier_dist = norm(b1 .- b0) + norm(b2 .- b1) + norm(b3 .- b2)
    
    return (cordal_dist + bezier_dist) * 0.5
end
"""
    estimate_catmullrom_arclength2d(p0, p1, p2, p3)

Uses 5-point Gauss-Legandre integration.
"""
function estimate_catmullrom_arclength2d(p0::T, p1::T, p2::T, p3::T) where T
    b0, b1, b2, b3 = catmullrom_as_bezier(p0, p1, p2, p3)
    c0 = complex(b0[1], b0[2])
    c1 = complex(b1[1], b1[2])
    c2 = complex(b2[1], b2[2])
    c3 = complex(b3[1], b3[2])
    
    return estimate_bezier_arclength(c0, c1, c2, c3)
end

#=
   refs
   https://qroph.github.io/2018/07/30/smooth-paths-using-catmull-rom-splines.html

   https://gist.github.com/njvack/6925609
    function catmullRom2bezier( points ) {
      // alert(points)
      var crp = points.split(/[,\s]/);

      var d = "";
      for (var i = 0, iLen = crp.length; iLen - 2 > i; i+=2) {
        var p = [];
        if ( 0 == i ) {
          p.push( {x: parseFloat(crp[ i ]), y: parseFloat(crp[ i + 1 ])} );
          p.push( {x: parseFloat(crp[ i ]), y: parseFloat(crp[ i + 1 ])} );
          p.push( {x: parseFloat(crp[ i + 2 ]), y: parseFloat(crp[ i + 3 ])} );
          p.push( {x: parseFloat(crp[ i + 4 ]), y: parseFloat(crp[ i + 5 ])} );
        } else if ( iLen - 4 == i ) {
          p.push( {x: parseFloat(crp[ i - 2 ]), y: parseFloat(crp[ i - 1 ])} );
          p.push( {x: parseFloat(crp[ i ]), y: parseFloat(crp[ i + 1 ])} );
          p.push( {x: parseFloat(crp[ i + 2 ]), y: parseFloat(crp[ i + 3 ])} );
          p.push( {x: parseFloat(crp[ i + 2 ]), y: parseFloat(crp[ i + 3 ])} );
        } else {
          p.push( {x: parseFloat(crp[ i - 2 ]), y: parseFloat(crp[ i - 1 ])} );
          p.push( {x: parseFloat(crp[ i ]), y: parseFloat(crp[ i + 1 ])} );
          p.push( {x: parseFloat(crp[ i + 2 ]), y: parseFloat(crp[ i + 3 ])} );
          p.push( {x: parseFloat(crp[ i + 4 ]), y: parseFloat(crp[ i + 5 ])} );
        }

        // Catmull-Rom to Cubic Bezier conversion matrix 
        //    0       1       0       0
        //  -1/6      1      1/6      0
        //    0      1/6      1     -1/6
        //    0       0       1       0

        // Cubic Bezier to Catmull-Rom conversion matrix inv(above)
        //   6       -6       0       1
        //   1        0       0       0
        //   0        0       0       1
        //   1        0      -6       6

        var bp = [];
        bp.push( { x: p[1].x,  y: p[1].y } );
        bp.push( { x: ((-p[0].x + 6*p[1].x + p[2].x) / 6), y: ((-p[0].y + 6*p[1].y + p[2].y) / 6)} );
        bp.push( { x: ((p[1].x + 6*p[2].x - p[3].x) / 6),  y: ((p[1].y + 6*p[2].y - p[3].y) / 6) } );
        bp.push( { x: p[2].x,  y: p[2].y } );

        d += "C" + bp[1].x + "," + bp[1].y + " " + bp[2].x + "," + bp[2].y + " " + bp[3].x + "," + bp[3].y + " ";
      }

      return d;
    }
   https://stackoverflow.com/questions/30748316/catmull-rom-interpolation-on-svg-paths
    For a curve segment defined by point P0, P1, P2 and P3 and knot sequence t0, t1, t2, t3, 
    the centripetal Catmull-Rom spline (defined between point P1 and P2) can be computed by
    the recursive formula provided in https://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline. Therefore, I will not elaborate here.

    To convert it to cubic Bezier curve, you need to compute the first derivative at P1 and P2 as

    M1 = (t2-t1)*(c1*(P1-P0)/(t1-t0) + c2*(P2-P1)/(t2-t1))
    M2 = (t2-t1)*(d1*(P2-P1)/(t2-t1) + d2*(P3-P2)/(t3-t2))
    Where

     c1 = (t2-t1)/(t2-t0),
     c2 = (t1-t0)/(t2-t0),
     d1 = (t3-t2)/(t3-t1),
     d2 = (t2-t1)/(t3-t1)

    Then you can convert it to cubic Bezier curve with 4 control points: Q0, Q1, Q2 and Q3:

    Q0 = P1
    Q1 = P1 + M1/3
    Q2 = P2 - M2/3
    Q3 = P2

=#
function catmullrom_as_bezier(p0::T, p1::T, p2::T, p3::T) where T
    b0 = p1
    b3 = p2
    
    d1 = norm(p1 .- p0); sqrtd1 = sqrt(d1)
    d2 = norm(p2 .- p1); sqrtd2 = sqrt(d2)
    d3 = norm(p3 .- p2); sqrtd3 = sqrt(d3)
    
    b1n = (d1 .* p2) - (d2 .* p0) + ((2*d1 + 3*sqrtd1*sqrtd2 + d2) .* p1)
    b1d = 3 * sqrtd1 * (sqrtd1 + sqrtd2)
    b1 = b1n ./ b1d
    b2n = (d3 .* p1) - (d2 .* p3) + ((2*d3 + 3*sqrtd3*sqrtd2 + d2) .* p2)
    b2d = 3 * sqrtd3 * (sqrtd3 + sqrtd2)
    b2 = b2n ./ b2d
    
    return b0, b1, b2, b3
end

# from https://github.com/fonttools/fonttools/blob/master/Lib/fontTools/misc/bezierTools.py

function estimate_bezier_arclength(c0::T, c1::T, c2::T, c3::T) where {T<:Complex}
   v0 = abs(c1-c0) * 0.15
   v1 = abs(-0.558983582205757*c0 + 0.325650248872424*c1 + 0.208983582205757*c2 + 0.024349751127576*c3)
   v2 = abs(c3-c0+c2-c1) * 0.26666666666666666
   v3 = abs(-0.024349751127576*c0 - 0.208983582205757*c1 - 0.325650248872424*c2 + 0.558983582205757*c3)
   v4 = abs(c3-c2) * 0.15
   return v0 + v1 + v2 + v3 + v4
end
