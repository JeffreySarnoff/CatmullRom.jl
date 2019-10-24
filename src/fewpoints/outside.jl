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

Base.convert(::Type{Array{T1,1}}, x::Array{T2,1}) where {N,T1,T2} = T1.(x)
Base.convert(::Type{NTuple{N,T1}}, x::NTuple{N,T2}) where {N,T1,T2} = T1.(x)
Base.convert(::Type{Array{NTuple{N,T},1}}, x::Array{T,1}) where {N,T} = map(a->(a...,), x)
Base.convert(::Type{Array{T,1}}, x::Array{NTuple{N,T},1}) where {N,T} = map(a->[a...,], x)
Base.convert(::Type{Array{T1,1}}, x::Array{NTuple{N,T2},1}) where {N,T1,T2} = convert(Array{T1,1}, convert(NTuple{N,T1},x))
Base.convert(::Type{Array{NTuple{N,T1},1}}, x::Array{T2,1}) where {N,T1,T2} = convert(NTuple{N,T1}, convert(Array{T1},x))

