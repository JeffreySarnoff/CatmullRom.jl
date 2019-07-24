using CatmullRom, Plots

b = 1/2;
a = 2*b;
fx(t) = a * cospi(t)^2;
fy(t) = (a^2/b) * cospi(t)^3 * sinpi(t);

n_points  =  7;
n_between = 17; # interpolants between each adjacent pair of points

xs = [fx(t) for t=range( 0.0, stop=1.0, length=n_points)];
ys = [fy(t) for t=range( 0.0, stop=1.0, length=n_points)];
points = collect(zip(xs, ys))
close_seq!(points)

crpoints = catmullrom(points, n_between);
cxs, cys = crpoints

plot(xs, ys, linecolor=:lightgreen,  linewidth=3, size=(600,600), legend=nothing, xaxis=nothing, yaxis=nothing)
plot!(cxs, cys, linecolor=:black, linewidth=2, size=(600,600), legend=nothing, xaxis=nothing, yaxis=nothing)


```

Here we use the same parameters on an open section of the Piriform curve.

```julia

n_points  =  9;
n_between = 24; # interpolants between each adjacent pair of points

xs = [fx(t) for t=range( -0.25, stop=0.25, length=n_points)];
ys = [fy(t) for t=range( -0.25, stop=0.25, length=n_points)];
points = collect(zip(xs,ys))

plot(points, linecolor=:blue,  size=(600,600), legend=nothing, xaxis=nothing, yaxis=nothing)
# overlay in fitted curve
cxs, cys = catmullrom(points, n_between);
plot!(cxs, cys, linecolor=:black, size=(600,600), legend=nothing, xaxis=nothing, yaxis=nothing)
```
