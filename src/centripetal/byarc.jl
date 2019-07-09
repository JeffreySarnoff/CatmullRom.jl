function catmullrombyarc(points::Points; arcpoints_min=2, arcpoints_max=64)
    xpoints = extendbounds(points);
    n_xpoints = npoints(xpoints)
    pointsperarc = arclength_interpolants(xpoints, arcpoints_min=arcpoints_min, arcpoints_max=arcpoints_max)
    n_points = sum(pointsperarc) + npoints(points)
    n_coords = ncoords(points)
    T = coordtype(points)
    result = fill(T[], n_coords)
    xpoint_seqs = [xpoints[i:i+3] for i=1:(length(xpoints)-3)]
    for idx=1:(n_xpoints-3)
        fourpoints = xpoint_seqs[idx]
        n_between  = pointsperarc[idx]
        fitted = catmullrom_splines(fourpoints, n_between)
        for j=1:n_coords
            append!(result[j], fitted[j][1:end-1])
        end    
    end
    for j=1:n_coords
        push!(result[j], points[end][j])
    end        
    return result
end
