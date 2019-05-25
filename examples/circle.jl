"""
Approaching the Unit Circle

    x(t) = cospi(t),    y(t) = sinpi(t)
    eigvecs(x, y) covers a circle as t ↦ 0..2
    
# We need 4 points to determine a Catmull-Rom spline over the _middle_ 2 points. 
# circleis a closed curve, so each pair of adjacent points is _middle_ to the other twts.

given₁ = [ (x= 1.0, y= 0.0), (x= 0.0, y= 1.0), (x= -1.0, y=0 .0), (x= 0.0, y= -1.0) ]

xs = first.(given₁), ys = final.(given₁)

Is that too coarse a list? Are there too few points given to support our purpose?

given₂ = [ (x=  1.0, y=  0.0), ( √3/2,  1/2), ( √2/2,  √2/2), ( 1/2,  √3/2), 
           (x=  0.0, y=  1.0), ( -1/2, √3/2), (-√2/2,  √2/2), (-√3/2,  1/2), 
           (x= -1.0, y=  0.0), (-√3/2, -1/2), (-√2/2, -√2/2), (-1/2, -√3/2),
           (x=  0.0, y= -1.0), ( √3/2,  1/2), ( √2/2,  √2/2), ( 1/2,  √3/2) ]

xsc, = given₁
ys = given₁


How many intermediating points should be introduced between adjacent point pairs?
"""
