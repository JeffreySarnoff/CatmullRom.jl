using CatmullRom, Plots

b = 1/2;
a = 2*b;
fx(t) = a * cospi(t)^2;
fy(t) = (a^2/b) * cospi(t)^3 * sinpi(t);

n_points  =  7;
n_between = 17; # interpolants between each adjacent pair of points

xs1 = [fx(t) for t=range( 0.0, stop=1.0, length=n_points)];
ys1 = [fy(t) for t=range( 0.0, stop=1.0, length=n_points)];
xs2 = reverse(xs1)
ys2 = reverse(ys1);

cxs1, cys1 = catmullrom(xs1, ys1, n_between);
cxs2, cys2 = catmullrom(xs2, ys2, n_between);

plot(xs2, ys2, linecolor=:lightgreen,  linewidth=3, size=(600,600), legend=nothing, xaxis=nothing, yaxis=nothing)
plot!(cxs2, cys2, linecolor=:black, linewidth=2, size=(600,600), legend=nothing, xaxis=nothing, yaxis=nothing)

#=
   The Piriform is a closed curve, so the first and last points are equal.
   Our two sequences of xs,ys have different starting points,
   so we need to pare off either the first or the last coordiantes
   before comparing the (sorted) sequences.
=#
sort(xs1[2:end]) == sort(xs2[2:end]) && sort(ys1[2:end]) == sort(ys2[2:end]) # true
  
#=
   As the two sequences of xs, ys contain the same coordiantes,
   it suffices to plot either one. Here we plot both to show
   that one exactly covers the other.
=#  
# plot(xs1, ys1, linecolor=:green, size=(600,600), legend=nothing, xaxis=nothing, yaxis=nothing)
plot(xs2, ys2, linecolor=:blue,  size=(600,600), legend=nothing, xaxis=nothing, yaxis=nothing)

# overlay in orange the curve fitted to the first sequence of xs,ys (xs1, ys1)
cxs1, cys1 = catmullrom(xs1, ys1, n_between);
plot!(cxs1, cys1, linecolor=:orange, size=(600,600), legend=nothing, xaxis=nothing, yaxis=nothing)

# overlay in black the curve fitted to the second sequence of xs,ys (xs2, ys2)
cxs2, cys2 = catmullrom(xs2, ys2, n_between);
plot!(cxs2, cys2, linecolor=:black, size=(600,600), legend=nothing, xaxis=nothing, yaxis=nothing)


```

There are two gaps between the orange and black curves, one at the upper right and the other at the lower left.
In the larger one, the orange curve is inside the black and in the smaller one the orange curve is outside the black.
For the larger gap, the black curve is more accurate.  For the smaller gap, the orange curve is more accurate.
The take-away message is that there will be some discrepancy when fitting closed curves with few given points.

Here we use the same parameters on an open section of the Piriform curve.

```julia
n_points  =  9;
n_between = 24; # interpolants between each adjacent pair of points

xs = [fx(t) for t=range( -0.25, stop=0.25, length=n_points)];
ys = [fy(t) for t=range( -0.25, stop=0.25, length=n_points)];

plot(xs, ys, linecolor=:blue,  size=(600,600), legend=nothing, xaxis=nothing, yaxis=nothing)
# overlay in fitted curve
cxs, cys = catmullrom(xs, ys, n_between);
plot!(cxs, cys, linecolor=:black, size=(600,600), legend=nothing, xaxis=nothing, yaxis=nothing)
```
