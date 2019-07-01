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

#=
zyxs = collect(zip(zs,ys,xs));
zs,ys,xs = catmullrom(collect(zip(zs,ys,xs)), 32, iterator=false);
zz,yy,xx = catmullrom(zyxs, 32, iterator=false);

plot(zz,yy,xx)
=#

CatmullRom.catmullrom(a::Base.Iterators.Flatten{Array{Array{T,1},1}}, n::Int; iterator=false) where {T} =
    catmullrom(collect(a), n, iterator=iterator)

CatmullRom.catmullrom(a::Base.Iterators.Zip{Tuple{Array{T,1},Array{T,1},Array{T,1}}}, n::Int; iterator=false) where {T} =
    catmullrom(collect(a), n, iterator=iterator)

zs,ys,xs = catmullrom(zip(zs,ys,xs), 32, iterator=false);
plot!(zs,ys,xs,linecolor=:navy)

```
