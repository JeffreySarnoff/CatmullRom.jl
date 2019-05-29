isclosed(firstpoint::V, lastpoint::V) where {T, V<:AbstractVector{T}} = 
    sum(abs.(lastpoint .- firstpoint)) < eps(T)^(5/8)

isclosed(firstpoint::NTuple{N,T}, lastpoint::NTuple{N,T}) where {N,T} = 
    sum(abs.(lastpoint .- firstpoint)) < eps(T)^(5/8)
    
isclosed(points) = isclosed(points[1], points[end])

augment_points(points) =
    isclosed(points[1], points[end]) ? augment_closed(points) : augment_open(points)

augment_closed(points) = [points[end-1], points..., points[2]]

augment_open(points) = [points[1], points..., points[end]]
