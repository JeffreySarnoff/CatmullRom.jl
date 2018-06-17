module Points

export Point, nd, setindex, Δpt², Δpt, 
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

#=
Point(x1::T) where {T} = Point((x,))
Point(x1::T, x2::T) where {T} = Point((x1, x2))
Point(x1::T, x2::T, x3::T) where {T} = Point((x1, x2, x3))
Point(x1::T, x2::T, x3::T, x4::T) where {T} = Point((x1, x2, x3, x4))
Point(x1::T, x2::T, x3::T, x4::T, x5::T) where {T} = Point((x1, x2, x3, x4, x5))
Point(x1::T, x2::T, x3::T, x4::T, x5::T, x6::T) where {T} = Point((x1, x2, x3, x4, x5, x6))
Point(x1::T, x2::T, x3::T, x4::T, x5::T, x6::T, x7::T) where {T} = Point((x1, x2, x3, x4, x5, x6, x7))
Point(x1::T, x2::T, x3::T, x4::T, x5::T, x6::T, x7::T, x8::T) where {T} = Point((x1, x2, x3, x4, x5, x6, x7, x8))
Point(x1::T, x2::T, x3::T, x4::T, x5::T, x6::T, x7::T, x8::T, x9::T) where {T} = Point((x1, x2, x3, x4, x5, x6, x7, x8, x9))
Point(x1::T, x2::T, x3::T, x4::T, x5::T, x6::T, x7::T, x8::T, x9::T, x10::T) where {T} = Point((x1, x2, x3, x4, x5, x6, x7, x8, x9, x10))
Point(x1::T, x2::T, x3::T, x4::T, x5::T, x6::T, x7::T, x8::T, x9::T, x10::T, x11::T) where {T} = Point((x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11))
Point(x1::T, x2::T, x3::T, x4::T, x5::T, x6::T, x7::T, x8::T, x9::T, x10::T, x11::T, x12::T) where {T} = Point((x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12))
Point(x1::T, x2::T, x3::T, x4::T, x5::T, x6::T, x7::T, x8::T, x9::T, x10::T, x11::T, x12::T, x13::T) where {T} = Point((x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13))
Point(x1::T, x2::T, x3::T, x4::T, x5::T, x6::T, x7::T, x8::T, x9::T, x10::T, x11::T, x12::T, x13::T, x14::T) where {T} = Point((x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14))
Point(x1::T, x2::T, x3::T, x4::T, x5::T, x6::T, x7::T, x8::T, x9::T, x10::T, x11::T, x12::T, x13::T, x14::T, x15::T) where {T} = Point((x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15))
Point(x1::T, x2::T, x3::T, x4::T, x5::T, x6::T, x7::T, x8::T, x9::T, x10::T, x11::T, x12::T, x13::T, x14::T, x15::T, x16::T) where {T} = Point((x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16))
=#

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

Δpt²(pt1::T) where {F, T<:Point{1,F}} = x1²(pt1)
Δpt²(pt1::T) where {F, T<:Point{2,F}} = x1²(pt1) + x2²(pt1) 
Δpt²(pt1::T) where {F, T<:Point{3,F}} = x1²(pt1) + x2²(pt1) + x3²(pt1) 
Δpt²(pt1::T) where {F, T<:Point{4,F}} = x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) 
Δpt²(pt1::T) where {F, T<:Point{5,F}} = x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) + x5²(pt1) 
Δpt²(pt1::T) where {F, T<:Point{6,F}} = x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) + x5²(pt1) + x6²(pt1) 
Δpt²(pt1::T) where {F, T<:Point{7,F}} = x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) + x5²(pt1) + x6²(pt1) + x7²(pt1) 
Δpt²(pt1::T) where {F, T<:Point{8,F}} = x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) + x5²(pt1) + x6²(pt1) + x7²(pt1) + x8²(pt1)
Δpt²(pt1::T) where {F, T<:Point{9,F}} = (x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) + x5²(pt1) + x6²(pt1) + x7²(pt1) + x8²(pt1)) + x9²(pt1)
Δpt²(pt1::T) where {F, T<:Point{10,F}} = (x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) + x5²(pt1) + x6²(pt1) + x7²(pt1) + x8²(pt1)) + (x9²(pt1) + x10²(pt1))
Δpt²(pt1::T) where {F, T<:Point{11,F}} = (x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) + x5²(pt1) + x6²(pt1) + x7²(pt1) + x8²(pt1)) + (x9²(pt1) + x10²(pt1) + x11²(pt1))
Δpt²(pt1::T) where {F, T<:Point{12,F}} = (x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) + x5²(pt1) + x6²(pt1) + x7²(pt1) + x8²(pt1)) + (x9²(pt1) + x10²(pt1) + x11²(pt1) + x12²(pt1))
Δpt²(pt1::T) where {F, T<:Point{13,F}} = (x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) + x5²(pt1) + x6²(pt1) + x7²(pt1) + x8²(pt1)) + (x9²(pt1) + x10²(pt1) + x11²(pt1) + x12²(pt1) + x13²(pt1))
Δpt²(pt1::T) where {F, T<:Point{14,F}} = (x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) + x5²(pt1) + x6²(pt1) + x7²(pt1) + x8²(pt1)) + (x9²(pt1) + x10²(pt1) + x11²(pt1) + x12²(pt1) + x13²(pt1) + x14²(pt1))
Δpt²(pt1::T) where {F, T<:Point{15,F}} = (x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) + x5²(pt1) + x6²(pt1) + x7²(pt1) + x8²(pt1)) + (x9²(pt1) + x10²(pt1) + x11²(pt1) + x12²(pt1) + x13²(pt1) + x14²(pt1) + x15²(pt1))
Δpt²(pt1::T) where {F, T<:Point{16,F}} = (x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) + x5²(pt1) + x6²(pt1) + x7²(pt1) + x8²(pt1)) + (x9²(pt1) + x10²(pt1) + x11²(pt1) + x12²(pt1) + x13²(pt1) + x14²(pt1) + x15²(pt1) + x16²(pt1))
Δpt²(pt1::T) where {F, T<:Point{17,F}} = (x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) + x5²(pt1) + x6²(pt1) + x7²(pt1) + x8²(pt1)) + (x9²(pt1) + x10²(pt1) + x11²(pt1) + x12²(pt1) + x13²(pt1) + x14²(pt1) + x15²(pt1) + x16²(pt1)) + (x17²(pt1))
Δpt²(pt1::T) where {F, T<:Point{18,F}} = (x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) + x5²(pt1) + x6²(pt1) + x7²(pt1) + x8²(pt1)) + (x9²(pt1) + x10²(pt1) + x11²(pt1) + x12²(pt1) + x13²(pt1) + x14²(pt1) + x15²(pt1) + x16²(pt1)) + (x17²(pt1) + x18²(pt1))
Δpt²(pt1::T) where {F, T<:Point{19,F}} = (x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) + x5²(pt1) + x6²(pt1) + x7²(pt1) + x8²(pt1)) + (x9²(pt1) + x10²(pt1) + x11²(pt1) + x12²(pt1) + x13²(pt1) + x14²(pt1) + x15²(pt1) + x16²(pt1)) + (x17²(pt1) + x18²(pt1) + x19²(pt1))
Δpt²(pt1::T) where {F, T<:Point{20,F}} = (x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) + x5²(pt1) + x6²(pt1) + x7²(pt1) + x8²(pt1)) + (x9²(pt1) + x10²(pt1) + x11²(pt1) + x12²(pt1) + x13²(pt1) + x14²(pt1) + x15²(pt1) + x16²(pt1)) + (x17²(pt1) + x18²(pt1) + x19²(pt1) + x20²(pt1))
Δpt²(pt1::T) where {F, T<:Point{21,F}} = (x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) + x5²(pt1) + x6²(pt1) + x7²(pt1) + x8²(pt1)) + (x9²(pt1) + x10²(pt1) + x11²(pt1) + x12²(pt1) + x13²(pt1) + x14²(pt1) + x15²(pt1) + x16²(pt1)) + (x17²(pt1) + x18²(pt1) + x19²(pt1) + x20²(pt1) + x21²(pt1))
Δpt²(pt1::T) where {F, T<:Point{22,F}} = (x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) + x5²(pt1) + x6²(pt1) + x7²(pt1) + x8²(pt1)) + (x9²(pt1) + x10²(pt1) + x11²(pt1) + x12²(pt1) + x13²(pt1) + x14²(pt1) + x15²(pt1) + x16²(pt1)) + (x17²(pt1) + x18²(pt1) + x19²(pt1) + x20²(pt1) + x21²(pt1) + x22²(pt1))
Δpt²(pt1::T) where {F, T<:Point{23,F}} = (x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) + x5²(pt1) + x6²(pt1) + x7²(pt1) + x8²(pt1)) + (x9²(pt1) + x10²(pt1) + x11²(pt1) + x12²(pt1) + x13²(pt1) + x14²(pt1) + x15²(pt1) + x16²(pt1)) + (x17²(pt1) + x18²(pt1) + x19²(pt1) + x20²(pt1) + x21²(pt1) + x22²(pt1) + x23²(pt1))
Δpt²(pt1::T) where {F, T<:Point{24,F}} = (x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) + x5²(pt1) + x6²(pt1) + x7²(pt1) + x8²(pt1)) + (x9²(pt1) + x10²(pt1) + x11²(pt1) + x12²(pt1) + x13²(pt1) + x14²(pt1) + x15²(pt1) + x16²(pt1)) + (x17²(pt1) + x18²(pt1) + x19²(pt1) + x20²(pt1) + x21²(pt1) + x22²(pt1) + x23²(pt1) + x24²(pt1))
Δpt²(pt1::T) where {F, T<:Point{25,F}} = (x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) + x5²(pt1) + x6²(pt1) + x7²(pt1) + x8²(pt1)) + (x9²(pt1) + x10²(pt1) + x11²(pt1) + x12²(pt1) + x13²(pt1) + x14²(pt1) + x15²(pt1) + x16²(pt1)) + (x17²(pt1) + x18²(pt1) + x19²(pt1) + x20²(pt1) + x21²(pt1) + x22²(pt1) + x23²(pt1) + x24²(pt1)) + x25²(pt1)
Δpt²(pt1::T) where {F, T<:Point{26,F}} = (x1²(pt1) + x2²(pt1) + x3²(pt1) + x4²(pt1) + x5²(pt1) + x6²(pt1) + x7²(pt1) + x8²(pt1)) + (x9²(pt1) + x10²(pt1) + x11²(pt1) + x12²(pt1) + x13²(pt1) + x14²(pt1) + x15²(pt1) + x16²(pt1)) + (x17²(pt1) + x18²(pt1) + x19²(pt1) + x20²(pt1) + x21²(pt1) + x22²(pt1) + x23²(pt1) + x24²(pt1)) + (x25²(pt1) + x26²(pt1))

const dpoint2 = Δpt²

# point, distance from the origin


Δpt(pt1::T) where {F, T<:Point{1,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{2,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{3,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{4,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{5,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{6,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{7,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{8,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{9,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{10,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{11,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{12,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{13,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{14,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{15,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{16,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{17,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{18,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{19,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{20,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{21,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{22,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{23,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{24,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{25,F}} = sqrt(Δpt²(pt1))
Δpt(pt1::T) where {F, T<:Point{26,F}} = sqrt(Δpt²(pt1))

const dpoint = Δpt


# separation for coordinate axes, oriented

x1(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x1(nt2) - x1(nt1)
x2(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x2(nt2) - x2(nt1)
x3(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x3(nt2) - x3(nt1)
x4(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x4(nt2) - x4(nt1)
x5(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x5(nt2) - x5(nt1)
x6(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x6(nt2) - x6(nt1)
x7(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x7(nt2) - x7(nt1)
x8(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x8(nt2) - x8(nt1)
x9(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x9(nt2) - x9(nt1)
x10(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x10(nt2) - x10(nt1)
x11(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x11(nt2) - x11(nt1)
x12(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x12(nt2) - x12(nt1)
x13(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x13(nt2) - x13(nt1)
x14(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x14(nt2) - x14(nt1)
x15(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x15(nt2) - x15(nt1)
x16(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x16(nt2) - x16(nt1)
x17(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x17(nt2) - x17(nt1)
x18(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x18(nt2) - x18(nt1)
x19(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x19(nt2) - x19(nt1)
x20(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x20(nt2) - x20(nt1)
x21(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x21(nt2) - x21(nt1)
x22(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x22(nt2) - x22(nt1)
x23(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x23(nt2) - x23(nt1)
x24(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x24(nt2) - x24(nt1)
x25(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x25(nt2) - x25(nt1)
x26(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = x26(nt2) - x26(nt1)

# squared separation for coordinate axes, unoriented

x1²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x1(nt1::T, nt2::T); d*d end
x2²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x2(nt1::T, nt2::T); d*d end
x3²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x3(nt1::T, nt2::T); d*d end
x4²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x4(nt1::T, nt2::T); d*d end
x5²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x5(nt1::T, nt2::T); d*d end
x6²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x6(nt1::T, nt2::T); d*d end
x7²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x7(nt1::T, nt2::T); d*d end
x8²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x8(nt1::T, nt2::T); d*d end
x9²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x9(nt1::T, nt2::T); d*d end
x10²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x10(nt1::T, nt2::T); d*d end
x11²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x11(nt1::T, nt2::T); d*d end
x12²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x12(nt1::T, nt2::T); d*d end
x13²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x13(nt1::T, nt2::T); d*d end
x14²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x14(nt1::T, nt2::T); d*d end
x15²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x15(nt1::T, nt2::T); d*d end
x16²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x16(nt1::T, nt2::T); d*d end
x17²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x17(nt1::T, nt2::T); d*d end
x18²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x18(nt1::T, nt2::T); d*d end
x19²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x19(nt1::T, nt2::T); d*d end
x20²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x20(nt1::T, nt2::T); d*d end
x21²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x21(nt1::T, nt2::T); d*d end
x22²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x22(nt1::T, nt2::T); d*d end
x23²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x23(nt1::T, nt2::T); d*d end
x24²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x24(nt1::T, nt2::T); d*d end
x25²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x25(nt1::T, nt2::T); d*d end
x26²(nt1::T, nt2::T) where {N, F, T<:Point{N,F}} = let d = x26(nt1::T, nt2::T); d*d end

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
