using CatmullRom, Plots
b = 1/2
a = 2*b
fx(t) = a * cospi(t)^2
fy(t) = (a^2/b) * cospi(t)^3 * sinpi(t)

xs = [fx(t) for t=0.0:0.02:1.0]
ys = [fy(t) for t=0.0:0.02:1.0]
xys = collect(zip(xs,ys))

if isapprox(xs[1], xs[end])
   xys = [xys[end-1], xys..., xys[2]]
end
cxs,cys = catmullrom(xys, 7)
plot(cxs,cys, size=(500,500))
  
