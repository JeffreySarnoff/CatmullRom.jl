
pa=[1,2];pb=[3,5];pc=[4,3];pd=[5,8];
pts=[pa,pb,pc,pd];

sqr(x) = x * x
distance(pa, pb) = sqrt(sum(sqr.(pb .- pa)))
dist(pa, pb) = sqrt(distance(pa, pb))

dists = [dist(pts[i],pts[i+1]) for i=1:length(pts)-1]
totaldist = sum(dists)
tidxs = [0.0, (cumsum(dists) ./ totaldist)...,]

