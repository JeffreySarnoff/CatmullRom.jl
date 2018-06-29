
pa=[1,2];pb=[3,5];pc=[4,3];pd=[5,8];

pts=[pa,pb,pc,pd];


sqr(x) = x * x

distance(pa, pb) = sqrt(sum(sqr.(pb .- pa)))

dist(pa, pb) = sqrt(distance(pa, pb))


dists = [dist(pts[i],pts[i+1]) for i=1:length(pts)-1]

dists = [dist(pts2D[i,:],pts2D[i+1,:]) for i=1:length(pts2D[:,1])-1]


totaldist = sum(dists)

tidxs = [0.0, (cumsum(dists) ./ totaldist)...,]


cheb(k,n)= (0+1)/2 + (1/2)*cospi(((2*(n+1-k)-1))/(2*n))


##############


pa=[1,2];pb=[3,5];pc=[4,3];pd=[5,8];

pts=[pa,pb,pc,pd];

distance1(pa, pb) = sqrt(sum(sqr.(pb .- pa)))

distance2(pa, pb) = sqrt(distance1(pa, pb))

dists(pts) = [distance2(pts[i],pts[i+1]) for i=1:(length(pts)-1)]

totaldist(pts) = sum(dists(pts))

tidxs(pts) = [0.0, (cumsum(dists(pts)) ./ totaldist(pts))...,]







ulia> pointss=[(0.0, 0.0),(1.0, 1.0),(2.0,1.0),(3.0, 0.0)]
4-element Array{Tuple{Float64,Float64},1}:
 (0.0, 0.0)
 (1.0, 1.0)
 (2.0, 1.0)
 (3.0, 0.0)

julia> interpolants=zero_chebroots_one(5)
7-element Array{Float64,1}:
 0.0                
 0.02447174185242318
 0.2061073738537635 
 0.5                
 0.7938926261462366 
 0.9755282581475768 
 1.0                

julia> tstchebroots = catmullrom((pointss...,), (interpolants...,))
7×2 Array{Float64,2}:
 1.0      1.0    
 1.02451  1.01187
 1.20856  1.07791
 1.50908  1.11051
 1.80333  1.06676
 1.97722  1.00924
 2.0      1.0    
