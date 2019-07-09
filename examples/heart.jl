
xt(t) = 16 * sin(t)^3
yt(t) = 13 * cos(t) - 5 * cos(2*t) - 2 * cos(3*t) - cos(4*t)

ts = collect(range(0.0,2.0*pi, length=15));
n_between = 28;

xs=xt.(ts); ys=yt.(ts);

crpts = catmullrom(xs,ys,n_between);
cxs,cys=crpts;

plot(xs, ys, linecolor=:lightgreen, linewidth=7, size=(500,500), xaxis=nothing, yaxis=nothing, legend=nothing)
plot!(cxs, cys, linecolor=:black, linewidth=2, size=(500,500), xaxis=nothing, yaxis=nothing, legend=nothing)

crpts = catmullrombyarc(collect(zip(xs,ys)), arcpoints_min=8, arcpoints_max=64);
cxs,cys=crpts;

plot(xs, ys, linecolor=:lightgreen, linewidth=7, size=(500,500), xaxis=nothing, yaxis=nothing, legend=nothing)
plot!(cxs, cys, linecolor=:black, linewidth=2, size=(500,500), xaxis=nothing, yaxis=nothing, legend=nothing)

#=

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
