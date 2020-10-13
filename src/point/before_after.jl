function pointbefore(firstpoints::P, scale=ReflectionScale) where P
    initialx = reflectback(first.(firstpoints[1:2])..., scale=scale)
    initialys = thiele4(firstpoints..., initialx)
    if isa(eltype(firstpoints), Tuple)
        initialpoint = (initialx, initialys[1:length(firstpoints[1])-1]...,)
    else
        initialpoint = [initialx, initialys[1:length(firstpoints[1])-1]...,]
    end    
    return convert(eltype(firstpoints), initialpoint)
end

function pointafter(lastpoints::P, scale=ReflectionScale) where P
    finalx = reflectforward(first.(lastpoints[end-1:end])..., scale=scale)
    finalys = thiele4(lastpoints..., finalx)
    finalpoint = (finalx, finalys...,)
    if isa(eltype(lastpoints), Tuple)
        finalpoint = (finalx, finalys[1:length(lastpoints[1])-1]...,)
    else
        finalpoint = [finalx, finalys[1:length(lastpoints[1])-1]...,]
    end    
    return convert(eltype(lastpoints), finalpoint)
end

allfinite(result) = all(Iterators.flatten(map(x->isfinite.(x),result)))
