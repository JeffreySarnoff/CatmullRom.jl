Approaching the Unit Circle

    x(t) = cospi(t);  y(t) = sinpi(t)
    (x(t), y(t)) circles as t ↦ 0..2

```julia

# We need 4 points to determine a Catmull-Rom spline over the _middle_ 2 points. 
# a circle is a closed curve, so the first and last points must be the same

given = [ (x=  1.0, y=  0.0), ( √3/2,  1/2), ( √2/2,  √2/2), ( 1/2,  √3/2), 
          (x=  0.0, y=  1.0), ( -1/2, √3/2), (-√2/2,  √2/2), (-√3/2,  1/2), 
          (x= -1.0, y=  0.0), (-√3/2, -1/2), (-√2/2, -√2/2), (-1/2, -√3/2),
          (x=  0.0, y= -1.0), ( √3/2,  1/2), ( √2/2,  √2/2), ( 1/2,  √3/2),
          (x=  1.0, y=  0.0)]

xs = first.(given₂); ys = last.(given₂)

# How many intermediating points should be introduced between adjacent point pairs?

nbetween = 8

plt = plot(title="MyPlot")

cxs,cys = catmullrom(Tuple.(given), nbetween);

plot!(first.(given), last.(given), size=(600,600))

plot!(cxs, cys, size=(600,600));

julia> plot!(cx,cy,size=(400,400));

julia> display(plt) #Display newly constructed plot

julia> cxs,cys = catmullrom(Tuple.(given₁), nbetween1^2);

julia> Plots.plot!(cxs,cys,size=(400,400));
=#

#=
 xs = first.(given₂); ys = last.(given₂)
 cxs,cys = catmullrom(Tuple.(given₂[4:10]), 17);
 cx=[sinpi(t) for t=1.5:0.02:2.0];
 cy=[cospi(t) for t=1.5:0.02:2.0];

plt = plot(title="MyPlot")
Plots.plot!(xs[5:9],ys[5:9],size=(400,400));
Plots.plot!(cxs,cys,size=(400,400));
display(plt) #Display newly constructed plot
Plots.plot!(cx,cy,size=(400,400));
display(plt) #Display newly constructed plot
=#
