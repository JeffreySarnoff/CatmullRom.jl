using CatmullRom, Plots

x(t) = cos(t)
y(t,m=4) = sin(t) * (sin(t/2)^m) # m=0..7

n_between = 24
given = collect(range(0.0, 2.0*pi, length=15));
xs = x.(given);
ys = y.(given, 5);

cmpts = catmullrom(xs,ys,n_between);
cxs,cys=cmpts;

plot(xs, ys, linecolor=:lightgreen, linewidth=9, size=(500,500), xaxis=nothing, yaxis=nothing, legend=nothing)
plot!(cxs, cys, linecolor=:black, linewidth=2, size=(500,500), xaxis=nothing, yaxis=nothing, legend=nothing)


crpts = catmullrombyarc(collect(zip(xs,ys)), arcpoints_min=6, arcpoints_max=16);
cxs,cys=crpts;

plot(xs, ys, linecolor=:lightgreen, linewidth=10, size=(500,500), xaxis=nothing, yaxis=nothing, legend=nothing)
plot!(cxs, cys, linecolor=:black, linewidth=2, size=(500,500), xaxis=nothing, yaxis=nothing, legend=nothing)
