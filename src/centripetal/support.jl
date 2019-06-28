isclosed(firstpoint::OnePoint, lastpoint::OnePoint) = firstpoint == lastpoint
isclosed(points::Points) = isclosed(first(points), last(points))

isopen(firstpoint::OnePoint, lastpoint::OnePoint) = firstpoint != lastpoint
isopen(points::Points) = isopen(first(points), last(points))


augment_points(points) =
    isclosed(points[1], points[end]) ? augment_closed(points) : augment_open(points)

augment_closed(points) = [points[end-1], points..., points[2]]

augment_open(points) = [points[1], points..., points[end]]
