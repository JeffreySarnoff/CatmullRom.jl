function catmullrombyarc(points::Points; arcpoints_min=2, arcpoints_max=64)
    xpoints = extendbounds(points);
    n_xpoints = npoints(xpoints)
    pointsperarc = arclength_interpolants(xpoints, arcpoints_min=arcpoints_min, arcpoints_max=arcpoints_max)
    n_points = sum(pointsperarc) + npoints(points)
    n_coords = ncoords(points)
    T = coordtype(points)
    result = Array{T,2}(undef, (n_points, n_coords))
    ptidx = 1
    for idx=1:(n_xpoints-3)
        fourpoints = xpoints[idx:idx+3]
        n_between  = pointsperarc[idx]
        fitted = catmullrom_splines(fourpoints, n_between)
        for j=1:n_coords
            result[ptidx:ptidx+n_between,j] = fitted[j][1:end-1][:]
        end    
        ptidx += n_between
    end
    return result
end
