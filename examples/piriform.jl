using CatmullRom, Plots

b = 1/2;
a = 2*b;
fx(t) = a * cospi(t)^2;
fy(t) = (a^2/b) * cospi(t)^3 * sinpi(t);

xs = [fx(t) for t=range(0.0, stop=1.0, length=16)];
ys = [fy(t) for t=range(0.0, stop=1.0, length=16)];

cxs, cys = catmullrom(xs, ys, 16);

plot(xs, ys, linecolor=:navy, size=(600,600), legend=nothing, xaxis=nothing, yaxis=nothing)
plot!(cxs,cys, linecolor=:green, size=(600,600), legend=nothing, xaxis=nothing, yaxis=nothing) 
