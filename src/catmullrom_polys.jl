
# compute Catmull-Rom cubic curve over [0, 1]
function catmullrom_cubic(x0::T, x1::T, x2::T, x3::T, dt0::T, dt1::T, dt2::T) where {T}
    # compute tangents when parameterized in [t1,t2]
    t1 = (x1 - x0) / dt0 - (x2 - x0) / (dt0 + dt1) + (x2 - x1) / dt1
    t2 = (x2 - x1) / dt1 - (x3 - x1) / (dt1 + dt2) + (x3 - x2) / dt2

    # rescale tangents for parametrization in [0,1]
    t1 *= dt1
    t2 *= dt1

    # return hermite cubic over [0,1]
    return hermite_cubic(x1, x2, t1, t2)
end

#=
   Compute coefficients for a cubic polynomial
     p(s) = c0 + c1*s + c2*s^2 + c3*s^3
   such that
     p(0) = x0, p(1) = x1
    and
     p'(0) = dx0, p'(1) = dx1.
=#
function hermite_cubic(x0::T, x1::T, dx0::T, dx1::T) where {T}
    c0 = x0
    c1 = dx0
    c2 = -3*x0 + 3*x1 - 2*dx0 - dx1
    c3 =  2*x0 - 2*x1 +   dx0 + dx1
    return Poly([c0, c1, c2, c3])
end

qrtrroot(x) = sqrt(sqrt(x))


function catmullrom_allpolys(points::Points{N,T}, deriv1::Bool=false, deriv2::Bool=false, integ1::Bool=false) where {N,T}
    polys = catmullrom_polys(points)

    d1polys = deriv1 ? polyder.(polys)    : nothing
    d2polys = deriv2 ? polyder.(d1polys)  : nothing
    i1polys = integ1 ? polyint.(polys)    : nothing

    result = (polys, d1polys, d2polys, i1polys)
    return result
end


#=
   given four x-ordinate sequenced ND points
   obtain N polys, parameterized over [0,1]
   interpolating from p1 to p2 inclusive
   one poly for each coordinate axis
=#
function catmullrom_polys(points::Points{N,T}) where {N,T}
    length(points) == 4 || throw(DomainError("exactly four points are required"))

    dt0, dt1, dt2 = prep_centripetal_catmullrom(points)
    pt0, pt1, pt2, pt3 = points

    dimen = length(pt0)
    polys = Vector{Poly{eltype(points[1])}}(undef, dimen)

    for i=1:dimen
        polys[i] = catmullrom_cubic(pt0[i], pt1[i], pt2[i], pt3[i], dt0, dt1, dt2)
    end

    return polys
end


#=
   determine the delta_traversal constants for the centripetal parameterization
      of the Catmull Rom cubic specified by four points (of increasing abcissae)
=#
function prep_centripetal_catmullrom(points::Points{N,T}) where {N,T}
    dt0 = qrtrroot(dot(points[1], points[2]))
    dt1 = qrtrroot(dot(points[2], points[3]))
    dt2 = qrtrroot(dot(points[3], points[4]))

    #safety check for repeated points
    if (dt1 < 1e-6) dt1 = 1.0 end
    if (dt0 < 1e-6) dt0 = dt1 end
    if (dt2 < 1e-6) dt2 = dt1 end

    return dt0, dt1, dt2
 end

