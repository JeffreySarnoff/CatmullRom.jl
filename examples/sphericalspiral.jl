http://mathworld.wolfram.com/SphericalSpiral.html


julia> ts=collect(range(-4pi,4pi, length=30));

julia> xs=x.(ts);

julia> ys=y.(ts);

julia> zs=z.(ts);

julia> plot(zs,ys,xs,linecolor=:navy)

julia> zyxs=collect(zip(zs,ys,xs));

julia> zz,yy,xx=catmullrom(zyxs, 32);

julia> plot(zz,yy,xx)
