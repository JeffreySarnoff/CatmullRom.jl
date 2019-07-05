x(t) = cos(t)
y(t,m=4) = sin(t) * (sin(t/2)^m) # m=0..7

n_between = 24
given = collect(range(0.0, 2.0*pi, length=15));
xs = x.(given);
ys = y.(given, 5);

plot(xs, ys)

cmpts = catmullrom(xs,ys,n_between);
cxs,cys=cmpts;
plot!(cxs,cys)
