
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


julia> tstchebroots01 = [into01(tstchebroots[:,1]) tstchebroots[:,2]]
7×2 Array{Float64,2}:
 0.0        1.0    
 0.0245142  1.01187
 0.208558   1.07791
 0.509085   1.11051
 0.803333   1.06676
 0.977221   1.00924
 1.0        1.0    
 
 
 julia> catmullrom((pointss...,), (tstchebroots[:,1]...,) )
7×2 Array{Float64,2}:
 1.0      1.0    
 1.02456  1.01189
 1.21106  1.07854
 1.51833  1.11021
 1.81256  1.06429
 1.9788   1.00861
 2.0      1.0    




julia> pointseq=[(0.0, 0.0),(1.0, 1.0),(2.0,1.0),(3.0, 0.0)]
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

julia> ccrpoints = catmullrom((pointseq...,), interpolants)
7×2 Array{Float64,2}:
 1.0      1.0    
 1.02451  1.01187
 1.20856  1.07791
 1.50908  1.11051
 1.80333  1.06676
 1.97722  1.00924
 2.0      1.0    


julia> ccrpullback = centripetal_pullback(ccrpoints)
14-element Array{Float64,1}:
 0.0                
 0.03611124334428015
 0.13505641482627284
 0.2614935746934738 
 0.3866032828406275 
 0.48277947012316846
 0.517589378111812  
 0.7482288197036723 
 0.7733555379702305 
 0.8326242506653972 
 0.8742727483720663 
 0.9225169699155896 
 0.9778327835664877 
 0.9999999999999998 


julia> into01(ans)
14-element Array{Float64,1}:
 0.0                 
 0.036111243344280154
 0.13505641482627287 
 0.26149357469347384 
 0.3866032828406276  
 0.4827794701231686  
 0.5175893781118122  
 0.7482288197036724  
 0.7733555379702307  
 0.8326242506653975  
 0.8742727483720665  
 0.9225169699155898  
 0.977832783566488   
 1.0 
 
 
 

