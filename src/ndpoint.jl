module Points

export Point, nd, setindex, Δpt², Δpt, dpoint2, dpoint,
    x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16,
    x17, x18, x19, x20, x21, x22, x23, x24, x25, x26,
    x1², x2², x3², x4², x5², x6², x7², x8²,  x9²,
    x10², x11², x12², x13², x14², x15², x16²,
    x17², x18², x19², x20², x21², x22², x23², x24², x25², x26²,
    x1, x2, x3, x4, x5, x6, x7, x8, x9,
    x10, x11, x12, x13, x14, x15, x16,
    x17, x18, x19, x20, x21, x22, x23, x24, x25, x26


# form points in 1D..26D coordinate space

# zs = zip(split(repeat("x",26),""), string.(collect(1:26)));
# coord_symbols = ([Symbol(string(i,j)) for (i,j) in zs]...,);
# (:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8, :x9, :x10, .., :x20, :x21, :x22, :x23, :x24, :x25, :x26)

struct Point{N,T}
    coords::NTuple{N,T}
end

Point(xs...) = Point((xs))

nd(x::Point{N,T}) where {T,N} = N
Base.eltype(x::Point{N,T}) where {T,N} = T

Base.lastindex(x::Point{N,T}) where {T,N} = N
Base.lastindex(::Type{Point{N,T}}) where {T,N} = N



Base.getindex(x::Point{N,T}, idx::I) where {T,N,I<:Union{Signed,Unsigned}} = x.coords[idx]
Base.getindex(x::Point{N,T}, idxs::R) where {T,N,R<:UnitRange} = x.coords[idxs]

function setindex(pt::Point{N,T}, value::T, idx::Signed) where {T,N}
    idx == 1 && return Point(value, pt[2:end]...,)
    idx == N && return Point(pt[1:end-1]..., value)
    return Point(pt[1:(idx-1)]...,value,pt[idx+1:end]...,)
end

function setindex(pt::Point{N,T}, values::NTuple{M,T}, idxs::R) where {T,N,M,R<:UnitRange}
    idxs.start == 1 && return Point(values..., pt[M+1:end]...,)
    idxs.stop == N && return Point(pt[1:end-M]..., values...,)
    return Point( pt[1:(idxs.start-1)]..., values..., pt[idxs.stop+1:end]...,)
end


x1(x::Point{N,T}) where {T,N} = x.coords[1]
x2(x::Point{N,T}) where {T,N} = x.coords[2]
x3(x::Point{N,T}) where {T,N} = x.coords[3]
x4(x::Point{N,T}) where {T,N} = x.coords[4]
x5(x::Point{N,T}) where {T,N} = x.coords[5]
x6(x::Point{N,T}) where {T,N} = x.coords[6]
x7(x::Point{N,T}) where {T,N} = x.coords[7]
x8(x::Point{N,T}) where {T,N} = x.coords[8]
x9(x::Point{N,T}) where {T,N} = x.coords[9]
x10(x::Point{N,T}) where {T,N} = x.coords[10]
x11(x::Point{N,T}) where {T,N} = x.coords[11]
x12(x::Point{N,T}) where {T,N} = x.coords[12]
x13(x::Point{N,T}) where {T,N} = x.coords[13]
x14(x::Point{N,T}) where {T,N} = x.coords[14]
x15(x::Point{N,T}) where {T,N} = x.coords[15]
x16(x::Point{N,T}) where {T,N} = x.coords[16]
x17(x::Point{N,T}) where {T,N} = x.coords[17]
x18(x::Point{N,T}) where {T,N} = x.coords[18]
x19(x::Point{N,T}) where {T,N} = x.coords[19]
x20(x::Point{N,T}) where {T,N} = x.coords[20]
x21(x::Point{N,T}) where {T,N} = x.coords[21]
x22(x::Point{N,T}) where {T,N} = x.coords[22]
x23(x::Point{N,T}) where {T,N} = x.coords[23]
x24(x::Point{N,T}) where {T,N} = x.coords[24]
x25(x::Point{N,T}) where {T,N} = x.coords[25]
x26(x::Point{N,T}) where {T,N} = x.coords[26]


# coordinate, squared distance from origin

x1²(x::Point{N,T}) where {T,N} = let a = x.coords[1]; a*a; end 
x2²(x::Point{N,T}) where {T,N} = let a = x.coords[2]; a*a; end
x3²(x::Point{N,T}) where {T,N} = let a = x.coords[3]; a*a; end
x4²(x::Point{N,T}) where {T,N} = let a = x.coords[4]; a*a; end
x5²(x::Point{N,T}) where {T,N} = let a = x.coords[5]; a*a; end
x6²(x::Point{N,T}) where {T,N} = let a = x.coords[6]; a*a; end
x7²(x::Point{N,T}) where {T,N} = let a = x.coords[7]; a*a; end
x8²(x::Point{N,T}) where {T,N} = let a = x.coords[8]; a*a; end
x9²(x::Point{N,T}) where {T,N} = let a = x.coords[9]; a*a; end
x10²(x::Point{N,T}) where {T,N} = let a = x.coords[10]; a*a; end
x11²(x::Point{N,T}) where {T,N} = let a = x.coords[11]; a*a; end
x12²(x::Point{N,T}) where {T,N} = let a = x.coords[12]; a*a; end
x13²(x::Point{N,T}) where {T,N} = let a = x.coords[13]; a*a; end
x14²(x::Point{N,T}) where {T,N} = let a = x.coords[14]; a*a; end
x15²(x::Point{N,T}) where {T,N} = let a = x.coords[15]; a*a; end
x16²(x::Point{N,T}) where {T,N} = let a = x.coords[16]; a*a; end
x17²(x::Point{N,T}) where {T,N} = let a = x.coords[17]; a*a; end
x18²(x::Point{N,T}) where {T,N} = let a = x.coords[18]; a*a; end
x19²(x::Point{N,T}) where {T,N} = let a = x.coords[19]; a*a; end
x20²(x::Point{N,T}) where {T,N} = let a = x.coords[20]; a*a; end
x21²(x::Point{N,T}) where {T,N} = let a = x.coords[21]; a*a; end
x22²(x::Point{N,T}) where {T,N} = let a = x.coords[22]; a*a; end
x23²(x::Point{N,T}) where {T,N} = let a = x.coords[23]; a*a; end
x24²(x::Point{N,T}) where {T,N} = let a = x.coords[24]; a*a; end
x25²(x::Point{N,T}) where {T,N} = let a = x.coords[25]; a*a; end
x26²(x::Point{N,T}) where {T,N} = let a = x.coords[26]; a*a; end

# point, squared distance from origin

Δpt²(pt::T) where {F, T<:Point{1,F}} = x1²(pt)
Δpt²(pt::T) where {F, T<:Point{2,F}} = x1²(pt) + x2²(pt) 
Δpt²(pt::T) where {F, T<:Point{3,F}} = x1²(pt) + x2²(pt) + x3²(pt) 
Δpt²(pt::T) where {F, T<:Point{4,F}} = x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) 
Δpt²(pt::T) where {F, T<:Point{5,F}} = x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) + x5²(pt) 
Δpt²(pt::T) where {F, T<:Point{6,F}} = x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) + x5²(pt) + x6²(pt) 
Δpt²(pt::T) where {F, T<:Point{7,F}} = x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) + x5²(pt) + x6²(pt) + x7²(pt) 
Δpt²(pt::T) where {F, T<:Point{8,F}} = x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) + x5²(pt) + x6²(pt) + x7²(pt) + x8²(pt)
Δpt²(pt::T) where {F, T<:Point{9,F}} = (x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) + x5²(pt) + x6²(pt) + x7²(pt) + x8²(pt)) + x9²(pt)
Δpt²(pt::T) where {F, T<:Point{10,F}} = (x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) + x5²(pt) + x6²(pt) + x7²(pt) + x8²(pt)) + (x9²(pt) + x10²(pt))
Δpt²(pt::T) where {F, T<:Point{11,F}} = (x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) + x5²(pt) + x6²(pt) + x7²(pt) + x8²(pt)) + (x9²(pt) + x10²(pt) + x11²(pt))
Δpt²(pt::T) where {F, T<:Point{12,F}} = (x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) + x5²(pt) + x6²(pt) + x7²(pt) + x8²(pt)) + (x9²(pt) + x10²(pt) + x11²(pt) + x12²(pt))
Δpt²(pt::T) where {F, T<:Point{13,F}} = (x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) + x5²(pt) + x6²(pt) + x7²(pt) + x8²(pt)) + (x9²(pt) + x10²(pt) + x11²(pt) + x12²(pt) + x13²(pt))
Δpt²(pt::T) where {F, T<:Point{14,F}} = (x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) + x5²(pt) + x6²(pt) + x7²(pt) + x8²(pt)) + (x9²(pt) + x10²(pt) + x11²(pt) + x12²(pt) + x13²(pt) + x14²(pt))
Δpt²(pt::T) where {F, T<:Point{15,F}} = (x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) + x5²(pt) + x6²(pt) + x7²(pt) + x8²(pt)) + (x9²(pt) + x10²(pt) + x11²(pt) + x12²(pt) + x13²(pt) + x14²(pt) + x15²(pt))
Δpt²(pt::T) where {F, T<:Point{16,F}} = (x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) + x5²(pt) + x6²(pt) + x7²(pt) + x8²(pt)) + (x9²(pt) + x10²(pt) + x11²(pt) + x12²(pt) + x13²(pt) + x14²(pt) + x15²(pt) + x16²(pt))
Δpt²(pt::T) where {F, T<:Point{17,F}} = (x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) + x5²(pt) + x6²(pt) + x7²(pt) + x8²(pt)) + (x9²(pt) + x10²(pt) + x11²(pt) + x12²(pt) + x13²(pt) + x14²(pt) + x15²(pt) + x16²(pt)) + (x17²(pt))
Δpt²(pt::T) where {F, T<:Point{18,F}} = (x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) + x5²(pt) + x6²(pt) + x7²(pt) + x8²(pt)) + (x9²(pt) + x10²(pt) + x11²(pt) + x12²(pt) + x13²(pt) + x14²(pt) + x15²(pt) + x16²(pt)) + (x17²(pt) + x18²(pt))
Δpt²(pt::T) where {F, T<:Point{19,F}} = (x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) + x5²(pt) + x6²(pt) + x7²(pt) + x8²(pt)) + (x9²(pt) + x10²(pt) + x11²(pt) + x12²(pt) + x13²(pt) + x14²(pt) + x15²(pt) + x16²(pt)) + (x17²(pt) + x18²(pt) + x19²(pt))
Δpt²(pt::T) where {F, T<:Point{20,F}} = (x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) + x5²(pt) + x6²(pt) + x7²(pt) + x8²(pt)) + (x9²(pt) + x10²(pt) + x11²(pt) + x12²(pt) + x13²(pt) + x14²(pt) + x15²(pt) + x16²(pt)) + (x17²(pt) + x18²(pt) + x19²(pt) + x20²(pt))
Δpt²(pt::T) where {F, T<:Point{21,F}} = (x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) + x5²(pt) + x6²(pt) + x7²(pt) + x8²(pt)) + (x9²(pt) + x10²(pt) + x11²(pt) + x12²(pt) + x13²(pt) + x14²(pt) + x15²(pt) + x16²(pt)) + (x17²(pt) + x18²(pt) + x19²(pt) + x20²(pt) + x21²(pt))
Δpt²(pt::T) where {F, T<:Point{22,F}} = (x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) + x5²(pt) + x6²(pt) + x7²(pt) + x8²(pt)) + (x9²(pt) + x10²(pt) + x11²(pt) + x12²(pt) + x13²(pt) + x14²(pt) + x15²(pt) + x16²(pt)) + (x17²(pt) + x18²(pt) + x19²(pt) + x20²(pt) + x21²(pt) + x22²(pt))
Δpt²(pt::T) where {F, T<:Point{23,F}} = (x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) + x5²(pt) + x6²(pt) + x7²(pt) + x8²(pt)) + (x9²(pt) + x10²(pt) + x11²(pt) + x12²(pt) + x13²(pt) + x14²(pt) + x15²(pt) + x16²(pt)) + (x17²(pt) + x18²(pt) + x19²(pt) + x20²(pt) + x21²(pt) + x22²(pt) + x23²(pt))
Δpt²(pt::T) where {F, T<:Point{24,F}} = (x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) + x5²(pt) + x6²(pt) + x7²(pt) + x8²(pt)) + (x9²(pt) + x10²(pt) + x11²(pt) + x12²(pt) + x13²(pt) + x14²(pt) + x15²(pt) + x16²(pt)) + (x17²(pt) + x18²(pt) + x19²(pt) + x20²(pt) + x21²(pt) + x22²(pt) + x23²(pt) + x24²(pt))
Δpt²(pt::T) where {F, T<:Point{25,F}} = (x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) + x5²(pt) + x6²(pt) + x7²(pt) + x8²(pt)) + (x9²(pt) + x10²(pt) + x11²(pt) + x12²(pt) + x13²(pt) + x14²(pt) + x15²(pt) + x16²(pt)) + (x17²(pt) + x18²(pt) + x19²(pt) + x20²(pt) + x21²(pt) + x22²(pt) + x23²(pt) + x24²(pt)) + x25²(pt)
Δpt²(pt::T) where {F, T<:Point{26,F}} = (x1²(pt) + x2²(pt) + x3²(pt) + x4²(pt) + x5²(pt) + x6²(pt) + x7²(pt) + x8²(pt)) + (x9²(pt) + x10²(pt) + x11²(pt) + x12²(pt) + x13²(pt) + x14²(pt) + x15²(pt) + x16²(pt)) + (x17²(pt) + x18²(pt) + x19²(pt) + x20²(pt) + x21²(pt) + x22²(pt) + x23²(pt) + x24²(pt)) + (x25²(pt) + x26²(pt))

const dpoipt2 = Δpt²


# point, distance from the origin


Δpt(pt::T) where {F, T<:Point{1,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{2,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{3,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{4,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{5,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{6,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{7,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{8,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{9,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{10,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{11,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{12,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{13,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{14,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{15,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{16,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{17,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{18,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{19,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{20,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{21,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{22,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{23,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{24,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{25,F}} = sqrt(Δpt²(pt))
Δpt(pt::T) where {F, T<:Point{26,F}} = sqrt(Δpt²(pt))


const dpoint = Δpt


# separation for coordinate axes, oriented

# separation for coordinate axes, oriented

x1(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x1(pt2) - x1(pt1)
x2(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x2(pt2) - x2(pt1)
x3(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x3(pt2) - x3(pt1)
x4(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x4(pt2) - x4(pt1)
x5(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x5(pt2) - x5(pt1)
x6(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x6(pt2) - x6(pt1)
x7(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x7(pt2) - x7(pt1)
x8(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x8(pt2) - x8(pt1)Δpt
x9(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x9(pt2) - x9(pt1)
x10(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x10(pt2) - x10(pt1)
x11(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x11(pt2) - x11(pt1)
x12(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x12(pt2) - x12(pt1)
x13(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x13(pt2) - x13(pt1)
x14(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x14(pt2) - x14(pt1)
x15(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x15(pt2) - x15(pt1)
x16(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x16(pt2) - x16(pt1)
x17(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x17(pt2) - x17(pt1)
x18(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x18(pt2) - x18(pt1)
x19(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x19(pt2) - x19(pt1)
x20(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x20(pt2) - x20(pt1)
x21(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x21(pt2) - x21(pt1)
x22(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x22(pt2) - x22(pt1)
x23(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x23(pt2) - x23(pt1)
x24(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x24(pt2) - x24(pt1)
x25(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x25(pt2) - x25(pt1)
x26(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = x26(pt2) - x26(pt1)

# squared separation for coordinate axes, unoriented

x1²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x1(pt1::T, pt2::T); d*d end
x2²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x2(pt1::T, pt2::T); d*d end
x3²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x3(pt1::T, pt2::T); d*d end
x4²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x4(pt1::T, pt2::T); d*d end
x5²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x5(pt1::T, pt2::T); d*d end
x6²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x6(pt1::T, pt2::T); d*d end
x7²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x7(pt1::T, pt2::T); d*d end
x8²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x8(pt1::T, pt2::T); d*d end
x9²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x9(pt1::T, pt2::T); d*d end
x10²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x10(pt1::T, pt2::T); d*d end
x11²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x11(pt1::T, pt2::T); d*d end
x12²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x12(pt1::T, pt2::T); d*d end
x13²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x13(pt1::T, pt2::T); d*d end
x14²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x14(pt1::T, pt2::T); d*d end
x15²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x15(pt1::T, pt2::T); d*d end
x16²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x16(pt1::T, pt2::T); d*d end
x17²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x17(pt1::T, pt2::T); d*d end
x18²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x18(pt1::T, pt2::T); d*d end
x19²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x19(pt1::T, pt2::T); d*d end
x20²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x20(pt1::T, pt2::T); d*d end
x21²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x21(pt1::T, pt2::T); d*d end
x22²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x22(pt1::T, pt2::T); d*d end
x23²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x23(pt1::T, pt2::T); d*d end
x24²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x24(pt1::T, pt2::T); d*d end
x25²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x25(pt1::T, pt2::T); d*d end
x26²(pt1::T, pt2::T) where {N, F, T<:Point{N,F}} = let d = x26(pt1::T, pt2::T); d*d end


# squared interpoint distance (norm2)

Δpt²(pt1::T, pt2::T) where {F, T<:Point{1,F}} = x1²(pt1, pt2)
Δpt²(pt1::T, pt2::T) where {F, T<:Point{2,F}} = x1²(pt1, pt2) + x2²(pt1, pt2) 
Δpt²(pt1::T, pt2::T) where {F, T<:Point{3,F}} = x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) 
Δpt²(pt1::T, pt2::T) where {F, T<:Point{4,F}} = x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) 
Δpt²(pt1::T, pt2::T) where {F, T<:Point{5,F}} = x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) + x5²(pt1, pt2) 
Δpt²(pt1::T, pt2::T) where {F, T<:Point{6,F}} = x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) + x5²(pt1, pt2) + x6²(pt1, pt2) 
Δpt²(pt1::T, pt2::T) where {F, T<:Point{7,F}} = x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) + x5²(pt1, pt2) + x6²(pt1, pt2) + x7²(pt1, pt2) 
Δpt²(pt1::T, pt2::T) where {F, T<:Point{8,F}} = x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) + x5²(pt1, pt2) + x6²(pt1, pt2) + x7²(pt1, pt2) + x8²(pt1, pt2)
Δpt²(pt1::T, pt2::T) where {F, T<:Point{9,F}} = (x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) + x5²(pt1, pt2) + x6²(pt1, pt2) + x7²(pt1, pt2) + x8²(pt1, pt2)) + x9²(pt1, pt2)
Δpt²(pt1::T, pt2::T) where {F, T<:Point{10,F}} = (x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) + x5²(pt1, pt2) + x6²(pt1, pt2) + x7²(pt1, pt2) + x8²(pt1, pt2)) + (x9²(pt1, pt2) + x10²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{11,F}} = (x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) + x5²(pt1, pt2) + x6²(pt1, pt2) + x7²(pt1, pt2) + x8²(pt1, pt2)) + (x9²(pt1, pt2) + x10²(pt1, pt2) + x11²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{12,F}} = (x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) + x5²(pt1, pt2) + x6²(pt1, pt2) + x7²(pt1, pt2) + x8²(pt1, pt2)) + (x9²(pt1, pt2) + x10²(pt1, pt2) + x11²(pt1, pt2) + x12²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{13,F}} = (x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) + x5²(pt1, pt2) + x6²(pt1, pt2) + x7²(pt1, pt2) + x8²(pt1, pt2)) + (x9²(pt1, pt2) + x10²(pt1, pt2) + x11²(pt1, pt2) + x12²(pt1, pt2) + x13²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{14,F}} = (x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) + x5²(pt1, pt2) + x6²(pt1, pt2) + x7²(pt1, pt2) + x8²(pt1, pt2)) + (x9²(pt1, pt2) + x10²(pt1, pt2) + x11²(pt1, pt2) + x12²(pt1, pt2) + x13²(pt1, pt2) + x14²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{15,F}} = (x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) + x5²(pt1, pt2) + x6²(pt1, pt2) + x7²(pt1, pt2) + x8²(pt1, pt2)) + (x9²(pt1, pt2) + x10²(pt1, pt2) + x11²(pt1, pt2) + x12²(pt1, pt2) + x13²(pt1, pt2) + x14²(pt1, pt2) + x15²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{16,F}} = (x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) + x5²(pt1, pt2) + x6²(pt1, pt2) + x7²(pt1, pt2) + x8²(pt1, pt2)) + (x9²(pt1, pt2) + x10²(pt1, pt2) + x11²(pt1, pt2) + x12²(pt1, pt2) + x13²(pt1, pt2) + x14²(pt1, pt2) + x15²(pt1, pt2) + x16²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{17,F}} = (x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) + x5²(pt1, pt2) + x6²(pt1, pt2) + x7²(pt1, pt2) + x8²(pt1, pt2)) + (x9²(pt1, pt2) + x10²(pt1, pt2) + x11²(pt1, pt2) + x12²(pt1, pt2) + x13²(pt1, pt2) + x14²(pt1, pt2) + x15²(pt1, pt2) + x16²(pt1, pt2)) + (x17²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{18,F}} = (x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) + x5²(pt1, pt2) + x6²(pt1, pt2) + x7²(pt1, pt2) + x8²(pt1, pt2)) + (x9²(pt1, pt2) + x10²(pt1, pt2) + x11²(pt1, pt2) + x12²(pt1, pt2) + x13²(pt1, pt2) + x14²(pt1, pt2) + x15²(pt1, pt2) + x16²(pt1, pt2)) + (x17²(pt1, pt2) + x18²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{19,F}} = (x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) + x5²(pt1, pt2) + x6²(pt1, pt2) + x7²(pt1, pt2) + x8²(pt1, pt2)) + (x9²(pt1, pt2) + x10²(pt1, pt2) + x11²(pt1, pt2) + x12²(pt1, pt2) + x13²(pt1, pt2) + x14²(pt1, pt2) + x15²(pt1, pt2) + x16²(pt1, pt2)) + (x17²(pt1, pt2) + x18²(pt1, pt2) + x19²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{20,F}} = (x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) + x5²(pt1, pt2) + x6²(pt1, pt2) + x7²(pt1, pt2) + x8²(pt1, pt2)) + (x9²(pt1, pt2) + x10²(pt1, pt2) + x11²(pt1, pt2) + x12²(pt1, pt2) + x13²(pt1, pt2) + x14²(pt1, pt2) + x15²(pt1, pt2) + x16²(pt1, pt2)) + (x17²(pt1, pt2) + x18²(pt1, pt2) + x19²(pt1, pt2) + x20²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{21,F}} = (x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) + x5²(pt1, pt2) + x6²(pt1, pt2) + x7²(pt1, pt2) + x8²(pt1, pt2)) + (x9²(pt1, pt2) + x10²(pt1, pt2) + x11²(pt1, pt2) + x12²(pt1, pt2) + x13²(pt1, pt2) + x14²(pt1, pt2) + x15²(pt1, pt2) + x16²(pt1, pt2)) + (x17²(pt1, pt2) + x18²(pt1, pt2) + x19²(pt1, pt2) + x20²(pt1, pt2) + x21²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{22,F}} = (x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) + x5²(pt1, pt2) + x6²(pt1, pt2) + x7²(pt1, pt2) + x8²(pt1, pt2)) + (x9²(pt1, pt2) + x10²(pt1, pt2) + x11²(pt1, pt2) + x12²(pt1, pt2) + x13²(pt1, pt2) + x14²(pt1, pt2) + x15²(pt1, pt2) + x16²(pt1, pt2)) + (x17²(pt1, pt2) + x18²(pt1, pt2) + x19²(pt1, pt2) + x20²(pt1, pt2) + x21²(pt1, pt2) + x22²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{23,F}} = (x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) + x5²(pt1, pt2) + x6²(pt1, pt2) + x7²(pt1, pt2) + x8²(pt1, pt2)) + (x9²(pt1, pt2) + x10²(pt1, pt2) + x11²(pt1, pt2) + x12²(pt1, pt2) + x13²(pt1, pt2) + x14²(pt1, pt2) + x15²(pt1, pt2) + x16²(pt1, pt2)) + (x17²(pt1, pt2) + x18²(pt1, pt2) + x19²(pt1, pt2) + x20²(pt1, pt2) + x21²(pt1, pt2) + x22²(pt1, pt2) + x23²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{24,F}} = (x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) + x5²(pt1, pt2) + x6²(pt1, pt2) + x7²(pt1, pt2) + x8²(pt1, pt2)) + (x9²(pt1, pt2) + x10²(pt1, pt2) + x11²(pt1, pt2) + x12²(pt1, pt2) + x13²(pt1, pt2) + x14²(pt1, pt2) + x15²(pt1, pt2) + x16²(pt1, pt2)) + (x17²(pt1, pt2) + x18²(pt1, pt2) + x19²(pt1, pt2) + x20²(pt1, pt2) + x21²(pt1, pt2) + x22²(pt1, pt2) + x23²(pt1, pt2) + x24²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{25,F}} = (x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) + x5²(pt1, pt2) + x6²(pt1, pt2) + x7²(pt1, pt2) + x8²(pt1, pt2)) + (x9²(pt1, pt2) + x10²(pt1, pt2) + x11²(pt1, pt2) + x12²(pt1, pt2) + x13²(pt1, pt2) + x14²(pt1, pt2) + x15²(pt1, pt2) + x16²(pt1, pt2)) + (x17²(pt1, pt2) + x18²(pt1, pt2) + x19²(pt1, pt2) + x20²(pt1, pt2) + x21²(pt1, pt2) + x22²(pt1, pt2) + x23²(pt1, pt2) + x24²(pt1, pt2)) + x25²(pt1, pt2)
Δpt²(pt1::T, pt2::T) where {F, T<:Point{26,F}} = (x1²(pt1, pt2) + x2²(pt1, pt2) + x3²(pt1, pt2) + x4²(pt1, pt2) + x5²(pt1, pt2) + x6²(pt1, pt2) + x7²(pt1, pt2) + x8²(pt1, pt2)) + (x9²(pt1, pt2) + x10²(pt1, pt2) + x11²(pt1, pt2) + x12²(pt1, pt2) + x13²(pt1, pt2) + x14²(pt1, pt2) + x15²(pt1, pt2) + x16²(pt1, pt2)) + (x17²(pt1, pt2) + x18²(pt1, pt2) + x19²(pt1, pt2) + x20²(pt1, pt2) + x21²(pt1, pt2) + x22²(pt1, pt2) + x23²(pt1, pt2) + x24²(pt1, pt2)) + (x25²(pt1, pt2) + x26²(pt1, pt2))


# interpoint distance

Δpt(pt1::T, pt2::T) where {F, T<:Point{1,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{2,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{3,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{4,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{5,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{6,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{7,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{8,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{9,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{10,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{11,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{12,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{13,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{14,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{15,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{16,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{17,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{18,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{19,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{20,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{21,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{22,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{23,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{24,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{25,F}} = sqrt(Δpt²(pt1, pt2))
Δpt(pt1::T, pt2::T) where {F, T<:Point{26,F}} = sqrt(Δpt²(pt1, pt2))

end # Points
