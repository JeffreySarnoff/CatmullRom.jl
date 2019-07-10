function thiele4(pt₁, pt₂, pt₃, pt₄, x)
    result = thiele4(pt₁[1], pt₂[1], pt₃[1], pt₄[1],
                     pt₁[2:end], pt₂[2:end], pt₃[2:end], pt₄[2:end], x)
    # check for degenerate case
    if !all(isfinite.(result))
        if x <= pt₃[1]
            result = thiele3(pt₁, pt₂, pt₃, x)
        else
            result = thiele3(pt₂, pt₃, pt₄, x)
        end    
    end
    return result
end

#=
function thiele4(pt₁x, pt₂x, pt₃x, pt₄x, pt₁y, pt₂y, pt₃y, pt₄y, x)
    t1 = @. pt₁y - pt₂y
    t2 = @. pt₂y - pt₃y
    t2 = @. inv(t2)
    t1 = @. inv(t1)
    t2 = @. (pt₂x - pt₃x) * t2
    t1 = @. (pt₁x - pt₂x) * t1
    t3 = @. -t2 + t1
    t4 = @. pt₃y - pt₄y
    t4 = @. inv(t4)
    t4 = @. -(pt₃x - pt₄x) * t4 + t2
    t4 = @. inv(t4)
    t3 = @. inv(t3)
    t3 = @. (pt₁x - pt₃x) * t3
    t4 = @. -(pt₂x - pt₄x) * t4 + pt₂y - pt₃y + t3
    t4 = @. inv(t4)
    t2 = @. -(pt₁x - pt₄x) * t4 + t1 - t2
    t2 = @. inv(t2)
    t2 = @. (pt₃x - x) * t2 - pt₁y + pt₂y + t3
    t2 = @. inv(t2)
    t1 = @. -(pt₂x - x) * t2 + t1
    t1 = @. inv(t1)
    return @. (-(pt₁x - x) * t1 + pt₁y)
end
=#

function thiele4(pt₁x, pt₂x, pt₃x, pt₄x, pt₁y, pt₂y, pt₃y, pt₄y, x)
    t1 = @. pt₁y - pt₂y
    t2 = @. pt₂y - pt₃y
    t2 = @. inv(t2)
    t1 = @. inv(t1)
    t2 =    (pt₂x .- pt₃x) .* t2
    t1 =    (pt₁x .- pt₂x) .* t1
    t3 = @. -t2 + t1
    t4 = @. pt₃y - pt₄y
    t4 = @. inv(t4)
    t4 =    -(pt₃x - pt₄x) .* t4 .+ t2
    t4 = @. inv(t4)
    t3 = @. inv(t3)
    t3 =    (pt₁x - pt₃x) .* t3
    t4 = @. -(pt₂x - pt₄x) * t4 + pt₂y - pt₃y + t3
    t4 = @. inv(t4)
    t2 = @. -(pt₁x - pt₄x) * t4 + t1 - t2
    t2 = @. inv(t2)
    t2 = @. (pt₃x - x) * t2 - pt₁y + pt₂y + t3
    t2 = @. inv(t2)
    t1 = @. -(pt₂x - x) * t2 + t1
    t1 = @. inv(t1)
    return @. (-(pt₁x - x) * t1 + pt₁y)
end
