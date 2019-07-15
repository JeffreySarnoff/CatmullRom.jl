using CatmullRom, Plots

x(t) = cos(t)
y(t,m=4) = sin(t) * (sin(t/2)^m) # m=0..7

n_between = 24
given = collect(range(0.0, length=15, stop=2.0*pi));
xs = x.(given);
ys = y.(given, 5);

points = collect(zip(xs,ys)
  
crpoints = catmullrom(points), n_between); #
cxs, cys = crpoints;                       # and
cxs, cys = catmullrom(points), n_between); # or


plot(xs, ys, linecolor=:lightgreen, linewidth=9, size=(500,500), xaxis=nothing, yaxis=nothing, legend=nothing)
plot!(cxs, cys, linecolor=:black, linewidth=2, size=(500,500), xaxis=nothing, yaxis=nothing, legend=nothing)


crpoints = catmullrom_byarc(collect(zip(xs,ys)), arcpoints_min=6, arcpoints_max=16);
cxs, cys = crpoints;

plot(xs, ys, linecolor=:lightgreen, linewidth=10, size=(500,500), xaxis=nothing, yaxis=nothing, legend=nothing)
plot!(cxs, cys, linecolor=:black, linewidth=2, size=(500,500), xaxis=nothing, yaxis=nothing, legend=nothing)
