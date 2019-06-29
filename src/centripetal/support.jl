isclosed(firstpoint::OnePoint, lastpoint::OnePoint) = firstpoint == lastpoint
isclosed(points::Points) = isclosed(first(points), last(points))

isopen(firstpoint::OnePoint, lastpoint::OnePoint) = firstpoint != lastpoint
isopen(points::Points) = isopen(first(points), last(points))

function myunzip(x)
    nd = length(x[1])
    npts = length(x)
    typ = eltype(eltype(x[1]))
    res = Array{typ, 2}(undef, (npts,nd))
    for k=1:npts
       for i=1:nd
          res[k,i] = x[k][i]
       end
    end
    return res
end
