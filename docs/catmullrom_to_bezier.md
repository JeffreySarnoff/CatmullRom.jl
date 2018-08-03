reference: 
    Glimpses of geometry and graphics
    Przemysław Koprowski

Let P1, P2, P3, P4 be the sequence of interpolation points for a Catmull-Rom segmental curve from P2 to P3.
Let Q1, Q2, Q3, Q4 be the sequence of control points for a Bezier segment that matches the Catmull-Rom segment.

Q1 = P2
Q2 = P1/(-6) + P2 + P3/6
Q3 = P2/6 + P3 + P4/(-6)
Q4 = P3
