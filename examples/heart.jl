#=

julia> given = collect(range(0.0,2.0*pi, length=12));

julia> n_between = 20;

julia> xs=x.(given);ys=y.(given);

julia> cmpts = catmullrom(xs,ys,n_between);

julia> plot(xs,ys)

julia> cxs,cys=cmpts;

julia> plot(cxs,cys)

julia> given = collect(range(0.0,2.0*pi, length=15));

julia> xs=x.(given);ys=y.(given);

julia> plot(xs,ys)

julia> cmpts = catmullrom(xs,ys,n_between);

julia> cxs,cys=cmpts;

julia> plot(cxs,cys)

julia> cmpts = catmullrom(xs,ys,32);

julia> cxs,cys=cmpts;

julia> plot(cxs,cys)
=#
