__precompile__()

module CentripetalCatmullRom

export CatmullRom, 
       PT1D, PT2D, PT3D, PT4D, coords,   # points in 1D, 2D, 3D, 4D [coordinates :x, :y, :z, :t]
       xcoord, ycoord, zcoord, tcoord,   # access point coordinates in 1D, 2D, 3D, 4D
       Δpoint, Δpoint2, dpoint, dpoint2, # dpoint, dpoint2 alias Δpoint, Δpoint2   
       polyval, polyder,                 # from Polynomials, specialized
       δpolyval, dpolyval                # dpolyval aliases δpolyval 

import Base: values                      # for NamedTuples

using Polynomials
import Polynomials: polyval, polyder


qrtrroot(x) = sqrt(sqrt(x))

# q.v. https://ideone.com/NoEbVM
#=
   Compute coefficients for a cubic polynomial
     p(s) = c0 + c1*s + c2*s^2 + c3*s^3
   such that
     p(0) = x0, p(1) = x1
    and
     p'(0) = dx0, p'(1) = dx1.
=#
function cubicpoly(x0::T, x1::T, dx0::T, dx1::T) where {T}
    c0 = x0
    c1 = dx0
    c2 = -3*x0 + 3*x1 - 2*dx0 - dx1
    c3 =  2*x0 - 2*x1 +   dx0 + dx1
    return Poly([c0, c1, c2, c3])
end

# compute coefficients for a nonuniform Catmull-Rom spline
function nonuniform_catmullrom(x0::T, x1::T, x2::T, x3::T, dt0::T, dt1::T, dt2::T) where {T}
    # compute tangents when parameterized in [t1,t2]
    t1 = (x1 - x0) / dt0 - (x2 - x0) / (dt0 + dt1) + (x2 - x1) / dt1
    t2 = (x2 - x1) / dt1 - (x3 - x1) / (dt1 + dt2) + (x3 - x2) / dt2
 
    # rescale tangents for parametrization in [0,1]
    t1 *= dt1
    t2 *= dt1
 
    return cubicpoly(x1, x2, t1, t2)
end

function prep_centripetal_catmullrom(p0::T, p1::T, p2::T, p3::T) where {N, F, T<:NTuple{N,F}}
    dt0 = qrtrroot(vecdot(p0, p1))
    dt1 = qrtrroot(vecdot(p1, p2))
    dt2 = qrtrroot(vecdot(p2, p3))
 
    #safety check for repeated points
    if (dt1 < 1e-4) dt1 = 1.0 end
    if (dt0 < 1e-4) dt0 = dt1 end
    if (dt2 < 1e-4) dt2 = dt1 end
 
    return dt0, dt1, dt2   
 end

function centripetal_catmullrom_polys(p0::T, p1::T, p2::T, p3::T) where {N, F, T<:NTuple{N,F}}
    dt0, dt1, dt2 = prep_centripetal_catmullrom(p0, p1, p2, p3)
 
    polys = Vector{Poly}(undef, N)
       
    for i=1:N
        polys[i] = nonuniform_catmullrom(p0[i], p1[i], p2[i], p3[i], dt0, dt1, dt2)
    end      
    
    return polys
end

function centripetal_catmullrom(points::Tuple{T,T,T,T}, interpolants::NTuple{M,F}) where {N, M, F, T<:NTuple{N,F}}
    polys = centripetal_catmullrom_polys(points...,)
    
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
       
end # module CentripetalCatmullRom
