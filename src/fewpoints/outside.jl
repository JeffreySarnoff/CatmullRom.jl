function pointbefore(firstpoints::P, scale=ReflectionScale) where P
    initialx = reflectback(first.(firstpoints[1:2])..., scale=scale)
    initialys = thiele4(firstpoints..., initialx)
    initialpoint = (initialx, initialys...,)
    return convert(eltype(firstpoints), initialpoint)
end

function pointafter(lastpoints::P, scale=ReflectionScale) where P
    finalx = reflectforward(first.(lastpoints[end-1:end])..., scale=scale)
    finalys = thiele4(lastpoints..., finalx)
    finalpoint = (finalx, finalys...,)
    return convert(eltype(lastpoints), finalpoint)
end

allfinite(result) = all(Iterators.flatten(map(x->isfinite.(x),result)))

Base.convert(::Type{Array{T,1}}, x::Array{NTuple{N,T},1} where {N,T} = map(a->[a...,], x)
Base.convert(::Type{Array{NTuple{N,T},1}}, x::Array{T,1} where {N,T} = map(a->(a...,), x)
