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
    return result
end

function quadraticBezierLength(p0, p1, p2)
    a = p0 .- 2 .* p1 .+ p2
    b = 2 .* p1 .- 2 .* p0
    A = 4 * sum(a .* a)
    B = 4 * sum(a .* b)
    C = sum(b .* b)

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
    return result
end

#=
robinrodricks commented on Apr 5, 2017
https://gist.github.com/tunght13488/6744e77c242cc7a94859

function quadraticBezierLength(p0, p1, p2) {
    var ax = p0.x - 2 * p1.x + p2.x;
    var ay = p0.y - 2 * p1.y + p2.y;
    var bx = 2 * p1.x - 2 * p0.x;
    var by = 2 * p1.y - 2 * p0.y;
    var A = 4 * (ax * ax + ay * ay);
    var B = 4 * (ax * bx + ay * by);
    var C = bx * bx + by * by;

    var Sabc = 2 * sqrt(A+B+C);
    var A_2 = sqrt(A);
    var A_32 = 2 * A * A_2;
    var C_2 = 2 * sqrt(C);
    var BA = B / A_2;

    return (A_32 * Sabc + A_2 * B * (Sabc - C_2) + (4 * C * A - B * B) * log((2 * A_2 + BA + Sabc) / (BA + C_2))) / (4 * A_32);
}
=#
