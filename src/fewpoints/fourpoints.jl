function thiele4(pt₁, pt₂, pt₃, pt₄, x)
    t1 = @. pt₁[2:end] - pt₂[2:end]
    t2 = @. pt₂[2:end] - pt₃[2:end]
    t2 = @. inv(t2)
    t1 = @. inv(t1)
    t2 = @. (pt₂[1] - pt₃[1]) * t2
    t1 = @. (pt₁[1] - pt₂[1]) * t1
    t3 = @. -t2 + t1
    t4 = @. pt₃[2:end] - pt₄[2:end]
    t4 = @. inv(t4)
    t4 = @. -(pt₃[1] - pt₄[1]) * t4 + t2
    t4 = @. inv(t4)
    t3 = @. inv(t3)
    t3 = @. (pt₁[1] - pt₃[1]) * t3
    t4 = @. -(pt₂[1] - pt₄[1]) * t4 + pt₂[2:end] - pt₃[2:end] + t3
    t4 = @. inv(t4)
    t2 = @. -(pt₁[1] - pt₄[1]) * t4 + t1 - t2
    t2 = @. inv(t2)
    t2 = @. (pt₃[1] - x) * t2 - pt₁[2:end] + pt₂[2:end] + t3
    t2 = @. inv(t2)
    t1 = @. -(pt₂[1] - x) * t2 + t1
    t1 = @. inv(t1)
    return @. (-(pt₁[1] - x) * t1 + pt₁[2:end])
end
