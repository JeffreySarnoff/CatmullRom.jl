It sounds like you want to keep the speed of the object at some constant value over the entire curve - knowing the arc-length won't help you do this. It will help you calculate at what time the object would reach its end-point if it were going at that speed, so it will be better than what you have now (the object will have the same average speed between all points), but the actual speed of the object will still vary as it moves around the curve.

A better solution would be to change our parametric parameter (the parameter that goes from 0 to 1, which I'll call s to avoid confusion with t = time) at a variable-rate ds/dt, which is determined by what speed you want the object to be moving at that point on the curve. So in other words, instead of changing s by 0.01 each frame, we might change it by 0.005 one frame, 0.02 the next, etc.

We do this by calculating the derivatives of x (dx/ds) and y (dy/ds) each frame, then setting

ds/dt = speed / sqrt( (dx/ds)2 + (dy/ds)2 )

That is, by taking the speed we want to go, and dividing by the speed we'd actually be going if we were changing s at a fixed increment.
Proof

We want the speed of our object to be constant; let's give that constant the name speed.

We learn in second-year calculus that, for parametric equations x(s) and y(s),

speed = sqrt( (dx/dt)2 + (dy/dt)2)

We also learn that

dx/dt = dx/ds * ds/dt    (chain rule)

Thus,

speed = sqrt( (dx/ds)2 (ds/dt)2 + (dy/ds)2 (ds/dt)2 )

Solving for ds/dt, we get the stated equation.
Calculating the derivatives

I've never worked with those particular splines, but I understand they just give x(s) and y(s) in terms of cubic-equations of s. Thus, we can find the derivative dx/ds easily: if

x(s) = a*s3 + b*s2 + c*s + e

then

dx/ds = 3a*s2 + 2b*s + c

(Same for dy/ds) Of course, you'll need to know the exact values of a, b, and c to do this. According this this page, those values are easy to find.
