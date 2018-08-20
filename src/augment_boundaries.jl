@inline function augmentends(::Type{Tuple}, endpoints::Function, points::Points{N,T}, closed) where {N,T}
    if !closed
        pre  = prepoint(endpoints, points[1:3]...,)
        post = postpoint(endpoints, points[end-2:end]...,)
    else
        pre  = points[end-2:end]
        post = points[1:3]
        #=
        pre   = (points[end] .+ points[1] .+ points[1])./3
        post  = (points[end] .+ points[end] .+ points[1])./3
        =#
    end

    return [(pre...,), points[1:end]..., (post...,)]
end

@inline function augmentends(::Type{Vector}, endpoints::Symbol, points::Points{N,T}, closed) where {N,T}
    if !closed
        pre  = prepoint(endpoints, points[1:3]...,)
        post = postpoint(endpoints, points[end-2:end]...,)
    else
        pre  = points[end-2:end]
        post = points[1:3]
    end

    return [(pre...,), points[1:end]..., (post...,)]
end



# julia> pre=CatmullRom.prepoint(xys[1:3]...,);

function prepoint(fn::F, point1::P, point2::P, point3::P) where {F<:Function, P<:OnePoint}
    dimen = length(point1)
    xpre = point1[1] - (point2[1] - point1[1])/8

    point = Vector{typeof(xpre)}(undef, dimen)
    point[1] = xpre
    p1 = point1[1:2]
    p2 = point2[1:2]
    p3 = point3[1:2]
    point[2] = fn(p1, p2, p3, xpre)

    for idx=3:dimen
        p1[2] = point1[idx]
        p2[2] = point2[idx]
        p3[2] = point3[idx]
        point[idx] = fn(p1, p2, p3, xpre)
    end

    return point
end


# julia> post=postpoint(xys[end-2:end]...,);

function postpoint(fn::F, point1::P, point2::P, point3::P) where {F<:Function, P<:OnePoint}
    dimen = length(point1)
    xpost = point3[1] + (point3[1] - point2[1])/8

    point = Vector{typeof(xpost)}(undef, dimen)
    point[1] = xpost
    p1 = point1[1:2]
    p2 = point2[1:2]
    p3 = point3[1:2]
    point[2] = fn(p1, p2, p3, xpost)

    for idx=3:dimen
        p1[2] = point1[idx]
        p2[2] = point2[idx]
        p3[2] = point3[idx]
        point[idx] = fn(p1, p2, p3, xpost)
    end

    return point
end
