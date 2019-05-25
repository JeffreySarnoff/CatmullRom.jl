firstone(xs) = xs[1]; lastone(xs) = xs[end]

"""
Approaching the Unit Circle

    x(t) = cospi(t);  y(t) = sinpi(t)
    (x(t), y(t)) circles as t ↦ 0..2
    
# We need 4 points to determine a Catmull-Rom spline over the _middle_ 2 points. 
# circleis a closed curve, so each pair of adjacent points is _middle_ to the other twts.

given₁ = [ (x= 1.0, y= 0.0), (x= 0.0, y= 1.0), (x= -1.0, y=0.0), (x= 0.0, y= -1.0) ]

xs = getfirst.(given₁); ys = getlast.(given₁)

Is that too coarse a list? Are there too few points given to support our purpose?

given₂ = [ (x=  1.0, y=  0.0), ( √3/2,  1/2), ( √2/2,  √2/2), ( 1/2,  √3/2), 
           (x=  0.0, y=  1.0), ( -1/2, √3/2), (-√2/2,  √2/2), (-√3/2,  1/2), 
           (x= -1.0, y=  0.0), (-√3/2, -1/2), (-√2/2, -√2/2), (-1/2, -√3/2),
           (x=  0.0, y= -1.0), ( √3/2,  1/2), ( √2/2,  √2/2), ( 1/2,  √3/2) ]

xs = getfirst.(given₁); ys = getlast.(given₁)
ys = getfirst.(given₂); getfirst.(


How many intermediating points should be introduced between adjacent point pairs? given₂?

