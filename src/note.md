
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


# mapping into centripetal curve 

sqr(x) = x * x
distance(pa, pb) = sqrt(sum(sqr.(pb .- pa)))
dist(pa, pb) = sqrt(distance(pa, pb))

# dists = [dist(pts[i],pts[i+1]) for i=1:length(pts)-1]

dists(pts) = [dist(pts[i,:],pts[i+1,:]) for i=1:length(pts[:,1])-1]

#totaldist = sum(dists)
centripetals(ptdists) = [0.0, (cumsum(ptdists) ./ sum(ptdists))...,]

centripetals(pts) = centripetals(dists(pts))


#########


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
 
 
 

julia> pointsa=[(points[i,:]...,) for i=1:7]
7-element Array{Tuple{Float64,Float64},1}:
 (1.0, 1.0)                              
 (1.0245142001653353, 1.0118687336747505)
 (1.208558366122857, 1.0779052011567294) 
 (1.5090845438775495, 1.1105137540145824)
 (1.8033334561990002, 1.0667591809605317)
 (1.9772207933291317, 1.0092375152298623)
 (2.0, 1.0)                              

julia> distances_centripetal(pointsa)
6-element Array{Float64,1}:
 0.16503406151018393
 0.44219089218288676
 0.5498091437072967 
 0.5454211768880193 
 0.42796546318842116
 0.15678318504323707

julia> ccrpullback = centripetal_pullback(pointsa)
7-element Array{Float64,1}:
 0.0                
 0.07215537708957281
 0.26548789450484556
 0.5058727322072556 
 0.7443390847339133 
 0.9314520303591936 
 1.0     
 
 
 julia> ccrpullback = centripetal_pullback(pointsa)
7-element Array{Float64,1}:
 0.0                
 0.07215537708957281
 0.26548789450484556
 0.5058727322072556 
 0.7443390847339133 
 0.9314520303591936 
 1.0                

julia> 

julia> points2=catmullrom((pointseq...,), ccrpullback)
7×2 Array{Float64,2}:
 1.0      1.0    
 1.07251  1.03291
 1.26925  1.0915 
 1.51506  1.11033
 1.75463  1.07873
 1.93577  1.02503
 2.0      1.0    

julia> diff(points)
┌ Warning: `diff(A::AbstractMatrix)` is deprecated, use `diff(A, dims=1)` instead.
│   caller = top-level scope at none:0
└ @ Core none:0
6×2 Array{Float64,2}:
 0.0245142   0.0118687 
 0.184044    0.0660365 
 0.300526    0.0326086 
 0.294249   -0.0437546 
 0.173887   -0.0575217 
 0.0227792  -0.00923752

julia> diff(diff(points))
┌ Warning: `diff(A::AbstractMatrix)` is deprecated, use `diff(A, dims=1)` instead.
│   caller = top-level scope at none:0
└ @ Core none:0
┌ Warning: `diff(A::AbstractMatrix)` is deprecated, use `diff(A, dims=1)` instead.
│   caller = top-level scope at none:0
└ @ Core none:0
5×2 Array{Float64,2}:
  0.15953      0.0541677
  0.116482    -0.0334279
 -0.00627727  -0.0763631
 -0.120362    -0.0137671
 -0.151108     0.0482842

julia> diff(points2)
┌ Warning: `diff(A::AbstractMatrix)` is deprecated, use `diff(A, dims=1)` instead.
│   caller = top-level scope at none:0
└ @ Core none:0
6×2 Array{Float64,2}:
 0.0725065   0.0329147
 0.196744    0.0585876
 0.245812    0.0188261
 0.239571   -0.0315946
 0.181141   -0.0537015
 0.0642257  -0.0250323

julia> diff(diff(points2))
┌ Warning: `diff(A::AbstractMatrix)` is deprecated, use `diff(A, dims=1)` instead.
│   caller = top-level scope at none:0
└ @ Core none:0
┌ Warning: `diff(A::AbstractMatrix)` is deprecated, use `diff(A, dims=1)` instead.
│   caller = top-level scope at none:0
└ @ Core none:0
5×2 Array{Float64,2}:
  0.124238     0.025673 
  0.0490683   -0.0397615
 -0.00624154  -0.0504207
 -0.0584299   -0.0221069
 -0.116915     0.0286692

 
 
 
 
 julia> xsinpoints=[k for k=0.0:(pi/8):pi]
9-element Array{Float64,1}:
 0.0                
 0.39269908169872414
 0.7853981633974483 
 1.1780972450961724 
 1.5707963267948966 
 1.9634954084936207 
 2.356194490192345  
 2.748893571891069  
 3.141592653589793  

julia> ysinpoints=[sin(k) for k=0.0:(pi/8):pi]
9-element Array{Float64,1}:
 0.0                   
 0.3826834323650898    
 0.7071067811865475    
 0.9238795325112867    
 1.0                   
 0.9238795325112867    
 0.7071067811865476    
 0.3826834323650899    
 1.2246467991473532e-16

julia> sinpoints =[ xsinpoints ysinpoints ]
9×2 Array{Float64,2}:
 0.0       0.0        
 0.392699  0.382683   
 0.785398  0.707107   
 1.1781    0.92388    
 1.5708    1.0        
 1.9635    0.92388    
 2.35619   0.707107   
 2.74889   0.382683   
 3.14159   1.22465e-16

julia> sinpointstup = [(sinpoints[i,:]...,) for i=1:length(sinpoints[:,1])]
9-element Array{Tuple{Float64,Float64},1}:
 (0.0, 0.0)                                 
 (0.39269908169872414, 0.3826834323650898)  
 (0.7853981633974483, 0.7071067811865475)   
 (1.1780972450961724, 0.9238795325112867)   
 (1.5707963267948966, 1.0)                  
 (1.9634954084936207, 0.9238795325112867)   
 (2.356194490192345, 0.7071067811865476)    
 (2.748893571891069, 0.3826834323650899)    
 (3.141592653589793, 1.2246467991473532e-16)


ulia> ccrpoints=catmullrom(sinpointstup,interpolants)
37×2 Array{Float64,2}:
 0.392699  0.382683
 0.459031  0.441703
 0.526422  0.500188
 0.593813  0.557055
 0.660145  0.61122 
 0.72436   0.661598
 0.785398  0.707107
 0.858687  0.757881
 0.927399  0.801222
 ⋮                 
 2.29239   0.749811
 2.35619   0.707107
 2.42339   0.659036
 2.48983   0.607407
 2.55557   0.553089
 2.62063   0.496949
 2.68506   0.439858
 2.74889   0.382683
