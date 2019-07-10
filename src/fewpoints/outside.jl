function pointbefore(firstpoints::Points, scale=ReflectionScale)
    initialx = reflectback(first.(firstpoints[1:2])..., scale=scale)
    initialys = thiele4(firstpoints..., initialx)
    initialpoint = (initialx, initialys...,)
    return convert(eltype(firstpoints), initialpoint)
end

function pointafter(lastpoints::Points, scaleReflectionScale)
    finalx = reflectforward(first.(lastpoints[end-1:end])..., scale=scale)
    finalys = thiele4(lastpoints..., finalx)
    finalpoint = (finalx, finalys...,)
    return convert(eltype(lastpoints), finalpoint)
end
