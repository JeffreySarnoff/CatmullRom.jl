


struct Cubic{T}
    coeffs::NTuple{4,T}
end
Cubic(a::T, b::T, c::T, d::T) where {T} = Cubic{T}((a,b,c,d))

#=
     Compute coefficients for a cubic polynomial
        p(s) = c0 + c1*s + c2*s^2 + c3*s^3
       such that
          p(0) = x0, p(1) = x1
        and
          p'(0) = t0, p'(1) = t1
=#
function HermiteCubic(x0::T, x1::T, t0::T, t1::T) where {T}
    c0 = x0
    c1 = t0
    c2 = -3 * x0 + 3 * x1 - 2 * t0 - t1
    c3 = 2 * x0 - 2 * x1 + t0 + t1
    return Cubic{T}((c0, c1, c2, c3))
end

# create alt version passing in Vector for results
function centripetal_catmullrom(p0::T, p1::T, p2::T, p3::T; abserr=1.0e-4) where {T}
    dt0 = root4(distancesquared(p0, p1))
    dt1 = root4(distancesquared(p1, p2))
    dt2 = root4(distancesquared(p2, p3))
	
    # safety check for repeated points
    small = eltype(T)(abserr)
    if dt1 < small
        dt1 = one(eltype(T))
    end
    if dt0 < small
        dt0 = dt1
    end
    if dt2 < small
        d2 = dt1
    end
	
    n = length(p0)
    result = Vector{Cubic}(undef, n)
    for i = 1:n
        result[i] = nonuniform_catmullrom(p0[i], p1[i], p2[i], p3[i], dt0, dt1, dt2)
    end
	
    return (result...,)	
end

#=
    nonuniform_catmullrom
    compute tangents when parameterized in [t1,t2]
    t1 = ( x1 - x0 ) / dt0 - ( x2 - x0 ) / ( dt0 + dt1 ) + ( x2 - x1 ) / dt1;
    t2 = ( x2 - x1 ) / dt1 - ( x3 - x1 ) / ( dt1 + dt2 ) + ( x3 - x2 ) / dt2;
    Maple simplify
    t1 = ((x2 - x1)*dt0^2 - dt1^2*(-x1 + x0)) / (dt0*(dt0 + dt1)*dt1)
    t2 = ((x3 - x2)*dt1^2 - dt2^2*(-x2 + x1)) / (dt1*(dt1 + dt2)*dt2)
    rescale tangents for parametrization in [0,1]
    t1 *= dt1
    t2 *= dt1
    Maple resimplify
    t1 = ((x2 - x1)*dt0^2 - dt1^2*(-x1 + x0)) / (dt0*(dt0 + dt1))
    t2 = ((x3 - x2)*dt1^2 - dt2^2*(-x2 + x1)) / ((dt1 + dt2)*dt2)
=#
    
function nonuniform_catmullrom(x0::T, x1::T, x2::T, x3::T, dt0::T, dt1::T, dt2::T) where {T}
    dt0sq = dt0 * dt0
    dt1sq = dt1 * dt1
    dt2sq = dt2 * dt2
    
    x2x1 = x2 - x1
    t1den = (dt0 + dt1) * dt0
    t2den = (dt1 + dt2) * dt2
    t1num = (x2x1 * dt0sq) - dt1sq * (x0 - x1)
    t2num = (x2x1 * dt2sq) - dt1sq * (x2 - x3)
    t1 = t1num / t1den
    t2 = t2num / t2den
    
    return HermiteCubic(x1, x2, t1, t2)
end

root4(x) = sqrt(sqrt(x))

distancesquared(a::T, b::T) where {T} = sum(subsq(a, b))

subsq(a::NTuple{N,T}, b::NTuple{N,T}) where {N,T} =
    (((a[i] - b[i])^2 for i=1:N)...,)

subsq(a::T, b::T) where {S,T<:AbstractVector{S}} =
    (((a[i] - b[i])^2 for i=1:min(length(a),length(b)))...,)

#


# Define a set of points for curve to go through
xcoords = Float64.([0, 2, 3, 4, 5, 6, 7])
ycoords = Float64.([1.5, 2, 1, 0.5, 1, 2, 3])
points = collect(zip(xcoords, ycoords))
npoints = length(points)
seqnums = collect(1:npoints)

nbetweenpoints = 4
nbetweens = collect(1:nbetweenpoints) ./ (nbetweenpoints + 1)

xs = [xcoords[1]]
ys = [ycoords[1]]
yofxs = [ycoords[1]]

for i=1:npoints-3
    seq2 = i+1
    seq3 = i+2
    dseq = seq3 - seq2
    seqbetweens = seq2 .+ nbetweens
    seqthrough  = (seq2, seqbetweens..., seq3)
    
    x2, y2 = points[i+1]
    x3, y3 = points[i+2]
    dx = x3 - x2
    dy = y3 - y2
    xinsides = x2 .+ (nbetweens .* dx)
    yinsides = y2 .+ (nbetweens .* dy)
    
    xpoly, ypoly = centripetal_catmullrom(points[i],points[i+1],points[i+2],points[i+3])
    xbetweens = (evalpoly(xinsides[idx], xpoly.coeffs) for idx=1:length(xinsides))
    ybetweens = (evalpoly(yinsides[idx], ypoly.coeffs) for idx=1:length(yinsides))
    yofxbetweens = (evalpoly(xbetweens[idx], ypoly.coeffs) for idx=1:length(xbetweens))
    
    append(xs, (xbetweens..., x3))
    append(ys, (ybetweens..., y3))
    append(yofxs, (yofxbetweens..., y3))
 end
 
 println(xs)
 println(ys)
 println(yofxs)
 
