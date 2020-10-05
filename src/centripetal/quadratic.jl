function catmullrom_to_approxarclength(p0, p1, p2, p3)
    bz0, bz1 = catmullrom_to_quadraticbeziers(p0, p1, p2, p3)
    result = quadraticBezierLength2d(bz0...)
    result += quadraticBezierLength2d(bz1...)
    return result
end

#=
Quadratic Approximation of Cubic Curves
NGHIA TRUONG, University of Utah
CEM YUKSEL, University of Utah
LARRY SEILER, Facebook Reality Lab
=#
function catmullrom_to_quadraticbeziers(p0, p1, p2, p3)
    p10 = p1 .- p0
    p21 = p2 .- p1
    p32 = p3 .- p2
    d1 = norm(p10); d1sqrt = sqrt(d1)
    d2 = norm(p21); d2sqrt = sqrt(d2)
    d3 = norm(p32); d3sqrt = sqrt(d3)

    q00 = p1
    q01 = p1 .+ ((d1 .* p21)+(d2 .* p10)) / (4*d1sqrt*(d1sqrt+d2sqrt))
    q11 = p2 .+ ((d3 .* (-1 .* p21)+(d2 .* (-1 .* p32))) / (4*d3sqrt*(d3sqrt+d2sqrt))
    q02 = q10 = (q01 .+ q11) ./ 2
    q12 = p2
    
    return (q00,q01,q02), (q10,q11,q12)
end

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
