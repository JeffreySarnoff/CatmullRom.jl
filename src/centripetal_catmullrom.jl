"""
     catmullrom(points, ninterpolants::Int; endpoints::Symbol=Thiele3, closed::Bool=false)
     catmullrom(points, interpolants; endpoints::Symbol=Thiele3, closed::Bool=false)

    `points` is a tuple/vector of points-as-tuples/vectors
    `interpolants` is a tuple/vector of values from 0.0 to 1.0 (inclusive)

all of the interpolants are applied to each segment

endpoints==Omit: interpolating from points[2] through points[end-1]
endpoints==Thiele3:  interpolating from points[1] through points[end]
"""
function catmullrom( points::Points{N,T},
                     interpolants::Seq01{M,S};
                     endpoints::Symbol=Thiele3,
                     closed::Bool=false
                   ) where {N,M,T<:Real,S<:AbstractFloat}

    npoints = length(points)
    npoints < 4 && throw(ErrorException("at least four points are required"))

    if endpoints != Omit
        TY = typeof(points[1]) <: Tuple ? Tuple : Vector
        points = augmentends(TY, endpointfn[endpoints], points, closed)
        npoints = length(points)
    end
    if npoints === 4
        catmullrom_4points(points, interpolants)
    else
        catmullrom_npoints(points, interpolants)
    end
end

function catmullrom_augmented( points::Points{N,T},
                               interpolants::Seq01{M,S};
                               endpoints::Symbol,
                               closed::Bool
                              ) where {N,M,T<:Real,S<:AbstractFloat}

     TY = typeof(points[1]) <: Tuple ? Tuple : Vector
     points = augmentends(TY, endpointfn[endpoints], points, closed)
     return catmullrom_npoints(points, interpolants)
end

function catmullrom(points::Points{N,T}, ninterpolants::Integer; endpoints::Symbol=Thiele3, closed::Bool=false) where {N,T}
    interpolants = uniform01(ninterpolants)
    return catmullrom(points, interpolants, endpoints=endpoints, closed=closed)
end


# ref https://ideone.com/NoEbVM


"""
    catmullrom_npoints(points, interpolants)
    `points` is a tuple of points-as-tuples
    `interpolants` is a tuple of values from 0.0 to 1.0 (inclusive)
interpolating points from points[2] through points[end-1] (inclusive)
"""
function catmullrom_npoints(pts::Points{N,T}, interpolants::Seq01{M,S}) where {N,T,M,S}
    npoints = length(pts)
    points_per_interpolation = length(interpolants)
    totalinterps = (npoints-4+1)*(points_per_interpolation - 1) + 1 # -1 for the shared end|1 point

    dimen = length(pts[1])
    TY = NTuple{dimen, typeof(pts[1][1])}
    points = Array{TY, 1}(undef, totalinterps)

    niters = npoints - 5

    points[1:points_per_interpolation, :] .= catmullrom_4points(pts[1:4], interpolants)

    idx₁ = 2; idx₂ = idx₁ + 3; sub₁ = 0; mul₁ = 1
    for k in 1:niters
        sub₂ = sub₁ - 1
        mul₂ = mul₁ + 1
        points[(mul₁ * points_per_interpolation + sub₁):(mul₂ * points_per_interpolation + sub₂), :] .=
            catmullrom_4points(pts[idx₁:idx₂], interpolants)
        idx₁ += 1; idx₂ += 1;
        sub₁, mul₁ = sub₂, mul₂
    end

    points[(end-points_per_interpolation+1):end, :] .= catmullrom_4points(pts[end-3:end], interpolants)

    return points
end


#=
   given four ND points in increasing abcissae order
   and k+2 interpolant values in 0..1 where 0.0 and 1.0 are the +2
   determine the k+2 interpolant determined ND points
   where the first interpolant point is the second ND point
   and the final interplant point is the third ND point
=#
function catmullrom_4points(pts::Points{N,T}, interpolants::Seq01{M,S}) where {N,T,M,S}
    polys = catmullrom_polys(pts)
    totalinterps = length(interpolants)
    dimen = length(pts[1])
    TY = typeof(pts[1][1])
    points = Array{TY, 2}(undef, (totalinterps,dimen))

    for col in 1:dimen
        points[1, col] = pts[2][col]
        points[end, col] = pts[3][col]
    end

    totalinterps -= 1

    for col in 1:dimen
        ply = polys[col]
        for row in 2:totalinterps
            value = interpolants[row]
            points[row, col] = polyval(ply, value)
        end
    end

    return [(points[i,:]...,) for i=1:size(points)[1]]
end

