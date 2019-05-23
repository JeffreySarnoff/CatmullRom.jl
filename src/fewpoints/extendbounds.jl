function extendbounds(extender::FUNCTION, 
                      firstpoints::T, finalpoints::T) where {T}
    return  firstpoint(extender, firstpoints[1:3]...,), 
            finalpoint(extender, finalpoints[1:3]...,)
end

function extendbounds(extender::FUNCTION,
                      points)    
    return firstpoint(extender, points[1:3]...,),
           finalpoint(extender, points[end-2:end]...,)
end

function extendbounds(firstextender::FUNCTION, finalextender::FUNCTION,
                      firstpoints, finalpoints)
    return  firstpoint(firstextender, firstpoints[1:3]...,), 
            finalpoint(finalextender, finalpoints[1:3]...,)
end

function extendbounds(firstextender::FUNCTION, finalextender::FUNCTION,
                      points)
    return firstpoint(firstextender, points[1:3]...,),
           finalpoint(finalextender, points[end-2:end]...,)
end

#=    
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
=#


function firstpoint(fn::F, point1::P, point2::P, point3::P) where {P, F<:Function}
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

function finalpoint(fn::F, point1::P, point2::P, point3::P) :: {P, F<:FUNCTION}
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
