Let P1, P2, P3, P4 be the sequence of interpolation points for a Catmull-Rom segmental curve from P2 to P3.
Let Q1, Q2, Q3, Q4 be the sequence of control points for a Bezier segment that matches the Catmull-Rom segment.


reference:
On the parameterization of Catmull-Rom curves - ‎Yuksel - Cited by 23
Parameterization and applications of Catmull–Rom … - ‎Yuksel - Cited by 43

The Centripetal Parmeterization (alpha = 1/2)
reparametrized into [0, 1]

Q1 = P2
Q2 = (d1\*P3 - d2\*P1 + ()\*P2) / denom
Q3 = (d3\*P2 - d2\*P3 + ()\*P3) / denom
Q4 = P3

denom = 3 * sqrt(d3) * (sqrt(d3) + sqrt(d2))

d1 = dist(P1,P2)
d2 = dist(P2,P3)    (the bounds of the Camull-Rom segment)
d3 = dist(P3,P4)



reference: 
    Glimpses of geometry and graphics
    Przemysław Koprowski


Q1 = P2
Q2 = P1/(-6) + P2 + P3/6
Q3 = P2/6 + P3 + P4/(-6)
Q4 = P3


