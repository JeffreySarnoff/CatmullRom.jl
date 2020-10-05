# https://gist.github.com/tunght13488/6744e77c242cc7a94859
function quadraticBezierLength2d(p0, p1, p2)
    ax = p0[1] - 2 * p1[1] + p2[1]
    ay = p0[2] - 2 * p1[2] + p2[2]
    bx = 2 * p1[1] - 2 * p0[1]
    by = 2 * p1[2] - 2 * p0[2]
    A = 4 * (ax * ax + ay * ay)
    B = 4 * (ax * bx + ay * by)
    C = bx * bx + by * by

    Sabc = 2 * sqrt(A+B+C)
    A2 = sqrt(A)
    A32 = 2 * A * A2
    C2 = 2 * sqrt(C)
    BA = B / A2

    r1 = A32 * Sabc
    r2 = A2 * B * (Sabc - C2)
    r3a = (4 * C * A - B * B)
    r3b = (2 * A2 + BA + Sabc) / (BA + C2)
    r3 = log(r3a * r3b)
    r4 = 4 * A32
    result = (r1 + r2 + r3) / r4
    return result # (A32 * Sabc + A2 * B * (Sabc - C2) + (4 * C * A - B * B) * log((2 * A2 + BA + Sabc) / (BA + C2))) / (4 * A32)
end
