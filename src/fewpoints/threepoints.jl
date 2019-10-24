#=
    Given three points (in their natural sequence) and a contained or nearby fourth abcissa,
    interpolate or extrapolate the fourth ordinate(s) using a polynomial or a rational approach.
    
        quadratic(pt₁, pt₂, pt₃, x)

        thiele3(pt₁, pt₂, pt₃, x)
=#

function quadratic(pt₁, pt₂, pt₃, x)
    result = quadratic(pt₁[1], pt₂[1], pt₃[1],
                       pt₁[2:end], pt₂[2:end], pt₃[2:end], x)
    # check for degenerate case
    if !allfinite.(result)
	result = linear(pt₁, pt₂, pt₃, x)
    end
	return result
end

function quadratic(pt₁x, pt₂x, pt₃x, pt₁y, pt₂y, pt₃y, x)
    t1 = @. pt₂x - pt₁x
    t2 = @. pt₁x - pt₃x
    t3 = @. pt₃x - pt₂x
    t4 = @. t2   * pt₂y
    t5 = @. pt₃x * pt₃x
    t6 = @. pt₂x * pt₂x
    t7 = @. pt₁x * pt₁x

    s = @. -(inv(t1) * inv(t2) * inv(t3))

    a = @. (pt₁y * t3 + pt₃y * t1 + t4) * x
    b = @. t5 * (pt₁y - pt₂y)
    c = @. t6 * (pt₃y - pt₁y)
    d = @. t7 * (pt₂y - pt₃y)
    q = @. t6 * (pt₁y * pt₃x - pt₃y * pt₁x)
    r = @. t4 *  pt₁x * pt₃x
    n = @. (-pt₁y * t5 + pt₃y * t7) * pt₂x
    aa = @. a - b - c - d
    bb = @. r - n - q

    result = @. (s * (aa * x + bb))
    return result # may contain !finite values
end


function thiele3(pt₁, pt₂, pt₃, x)
    result = thiele3(pt₁[1], pt₂[1], pt₃[1],
                     pt₁[2:end], pt₂[2:end], pt₃[2:end], x)
    # check for degenerate case
    if !allfinite.(result)
        result = quadratic(pt₁, pt₂, pt₃, x)
    end
    return result
end

function thiele3(pt₁x, pt₂x, pt₃x, pt₁y, pt₂y, pt₃y, x)
    t1 = @. pt₁y - pt₂y
    t2 = @. pt₂y - pt₃y
    t1 = @. inv(t1)
    t2 = @. inv(t2)
    t1 = @. (pt₁x - pt₂x) * t1
    t2 = @. -(pt₂x - pt₃x) * t2 + t1
    t2 = @. inv(t2)
    t2 = @. (pt₁x - pt₃x) * t2 - pt₁y + pt₂y
    t2 = @. inv(t2)
    t1 = @. -(pt₂x - x) * t2 + t1
    t1 = @. inv(t1)

    result = @. (-(pt₁x - x) * t1 + pt₁y)
    if !allfinite(result)
        result = (pt₁x, pt₂x, pt₃x, pt₁y, pt₂y, pt₃y, x)
    end

    return result
end
