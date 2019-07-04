"""
Approaching the Unit Circle

    x(t) = cospi(t);  y(t) = sinpi(t)
    (x(t), y(t)) circles as t ↦ 0..2
    
# We need 4 points to determine a Catmull-Rom spline over the _middle_ 2 points. 
# circleis a closed curve, so each pair of adjacent points is _middle_ to the other twts.

given₁ = [ (x= 1.0, y= 0.0), (x= 0.0, y= 1.0), (x= -1.0, y=0.0), (x= 0.0, y= -1.0) ]

xs = first.(given₁); ys = last.(given₁)

Is that too coarse a list? Are there too few points given to support our purpose?

given₂ = [ (x=  1.0, y=  0.0), ( √3/2,  1/2), ( √2/2,  √2/2), ( 1/2,  √3/2), 
           (x=  0.0, y=  1.0), ( -1/2, √3/2), (-√2/2,  √2/2), (-√3/2,  1/2), 
           (x= -1.0, y=  0.0), (-√3/2, -1/2), (-√2/2, -√2/2), (-1/2, -√3/2),
           (x=  0.0, y= -1.0), ( √3/2,  1/2), ( √2/2,  √2/2), ( 1/2,  √3/2) ]

xs = first.(given₂); ys = last.(given₂)


How many intermediating points should be introduced between adjacent point pairs? given₂?

nbetween = 3
=#

#=

julia> plt = plot(title="MyPlot")

julia> cxs,cys = catmullrom(Tuple.(given₁), 3);

julia> plot!(first.(given₁)[2:end-1],last.(given₁)[2:end-1],size=(400,400))

julia> plot!(cxs,cys,size=(400,400));

julia> plot!(cx,cy,size=(400,400));

julia> display(plt) #Display newly constructed plot

julia> cxs,cys = catmullrom(Tuple.(given₁), 17);

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
