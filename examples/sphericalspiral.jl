# http://mathworld.wolfram.com/SphericalSpiral.html

```julia
using CatmullRom, Plots

x(t) = cos(t) / sqrt(1 + t*t)
y(t) = sin(t) / sqrt(1 + t*t)
z(t) = - t/sqrt(1+t*t)
                        
ts = collect(range(-4pi,4pi, length=24));
xs = x.(ts); ys = y.(ts); zs = z.(ts);

plot(zs, ys, xs, linecolor=:darkred, linewidth=2, size=(600,600))

czs, cys, cxs = catmullrom(collect(zip(zs,ys,xs)), 36);
plot!(czs, cys, cxs, linecolor=:navy, linewidth=2, size=(600,600))

```
