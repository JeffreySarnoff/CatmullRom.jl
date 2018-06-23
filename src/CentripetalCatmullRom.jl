__precompile__()

module CentripetalCatmullRom

export catmullrom, into01

using Polynomials
import Polynomials: polyval, polyder

using LinearAlgebra: dot

dpolyval(ply::Poly, value::T) where {T} = polyval(polyder(ply), value)


"""
    into01((xs...,))
    into01([xs...,])

maps values into 0.0:1.0 (minimum(xs) --> 0.0, maximum(xs) --> 1.0)
"""
function into01(values::U) where {N,T, U<:Union{NTuple{N,T}, Vector{T}}}
    mn, mx = minimum(values), maximum(values)
    delta = mx - mn
    map(clamp01, (values .- mn) ./ delta)
end

@inline clamp01(x::T) where {T<:Real} = clamp(x, zero(T), one(T))
              
qrtrroot(x) = sqrt(sqrt(x))

# ref https://ideone.com/NoEbVM

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

# compute nonuniform Catmull-Rom spline over [0, 1]
function nonuniform_catmullrom(x0::T, x1::T, x2::T, x3::T, dt0::T, dt1::T, dt2::T) where {T}
    # compute tangents when parameterized in [t1,t2]
    t1 = (x1 - x0) / dt0 - (x2 - x0) / (dt0 + dt1) + (x2 - x1) / dt1
    t2 = (x2 - x1) / dt1 - (x3 - x1) / (dt1 + dt2) + (x3 - x2) / dt2
 
    # rescale tangents for parametrization in [0,1]
    t1 *= dt1
    t2 *= dt1
 
    # return hermite cubic over [0,1]
    return hermite_cubic(x1, x2, t1, t2)
end

function prep_catmullrom(p0::T, p1::T, p2::T, p3::T) where {N, F, T<:NTuple{N,F}}
    dt0 = qrtrroot(dot(p0, p1))
    dt1 = qrtrroot(dot(p1, p2))
    dt2 = qrtrroot(dot(p2, p3))
 
    #safety check for repeated points
    if (dt1 < 1e-4) dt1 = 1.0 end
    if (dt0 < 1e-4) dt0 = dt1 end
    if (dt2 < 1e-4) dt2 = dt1 end
 
    return dt0, dt1, dt2   
 end

#=
   given four x-ordinate sequenced ND points
   obtain N polys, parameterized over [0,1]
   interpolating from p1 to p2 inclusive
   one poly for each coordinate axis
=#
function catmullrom_polys(p0::T, p1::T, p2::T, p3::T) where {N, F, T<:NTuple{N,F}}
    dt0, dt1, dt2 = prep_catmullrom(p0, p1, p2, p3)
 
    polys = Vector{Poly}(undef, N)
       
    for i=1:N
        polys[i] = nonuniform_catmullrom(p0[i], p1[i], p2[i], p3[i], dt0, dt1, dt2)
    end      
    
    return polys
end


function catmullrom_points(pta::T, pt0::T, pt1::T, ptb::T, interpolants::NTuple{M,F}) where {N, M, F, T<:NTuple{N,F}}
    polys = catmullrom_polys(pta, pt0, pt1, ptb)
    
    points = Array{F, 2}(undef, (M,N))
    for col in 1:N
        ply = polys[col]
        for row in 1:M
            value = interpolants[row]
            0.0 <= value <= 1.0 || throw(DomainError("interpolant value ($value) should be in 0.0:1.0"))
            points[row, col] = polyval(ply, value)
        end
    end
    
    return points
end

#=
    interpolants
    1 2 3 4 5 6         1 2 3 4 5 6
    0         1         1         1
              1 2 3 4 5 6         1 2 3 4 5 6
    1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 6

    (ninterpolants-1)*(npointspans-1) + ninterpolants
=#
function catmullrom(points::U1, interpolants::U2) where {U1, U2}
    npoints = length(points)
    npoints < 4 && throw(ErrorException("at least four points are required"))
    point_windows = npoints - 3
    ninterp = length(interpolants)
    ninterp < 2 && throw(ErrorException("at least two interpolants [0,1] are required"))
    ndimens = length(points[1])
    eltyp = eltype(points[1])
    result = Array{eltyp, ndimens}(undef, (point_windows * ninterp, ndimens) )
    interpolants1 = view(interpolants,1:(ninterp-1))
    for i in 1:(point_windows-1)
         result[i:i+(ninterp-1),:] = catmullrom_points(points[i:i+3]..., interpolants1)
    end
    result[point_windows:(point_windows+ninterp), :] = catmullrom_points(points[end-3:end]..., interpolants)

    return result
end



end # module CentripetalCatmullRom
