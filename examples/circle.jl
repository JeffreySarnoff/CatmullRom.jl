Approaching the Unit Circle

    x(t) = cospi(t);  y(t) = sinpi(t)
    (x(t), y(t)) circles as t ↦ 0..2

```julia
using CatmullRom, Plots

fx(t) = cospi(t);
fy(t) = sinpi(t);

xs = [fx(t) for t=range(0.0, stop=2.0, length=17)];
ys = [fy(t) for t=range(0.0, stop=2.0, length=17)];

cxs,cys = catmullrom(xs, ys, 16);

plot(xs, ys, linecolor=:black, size=(600,600), legend=nothing, xaxis=nothing)
plot!(cxs,cys, linecolor=:blue, size=(600,600), legend=nothing, yaxis=nothing) 
```

