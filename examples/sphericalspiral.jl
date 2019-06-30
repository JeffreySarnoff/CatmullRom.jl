# http://mathworld.wolfram.com/SphericalSpiral.html

```julia
using CatmullRom, Plots

x(t) = cos(t) / sqrt(1 + t*t)
y(t) = sin(t) / sqrt(1 + t*t)
z(t) = - t/sqrt(1+t*t)
                        
ts = collect(range(-4pi,4pi, length=30));
xs = x.(ts);
ys = y.(ts);
zs = z.(ts);

plot(zs,ys,xs,linecolor=:navy)

zyxs = collect(zip(zs,ys,xs));
zz,yy,xx = catmullrom(zyxs, 32);

plot(zz,yy,xx)
```
