function reflectback(pt₁::T, pt₂::T; scale=1.0) where {T}
   @. pt₁ - ((pt₂ - pt₁) * scale)
end

function reflectforward(pt₁::T, pt₂::T; scale=1.0) where {T}
    return @. pt₂ - ((pt₁ - pt₂) * scale)
end
