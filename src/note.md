
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
