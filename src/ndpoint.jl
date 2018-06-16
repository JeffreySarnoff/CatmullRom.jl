# form points in 1D..26D coordinate space

# zs = zip(split(repeat("x",26),""), string.(collect(1:26)));
# coord_symbols = ([Symbol(string(i,j)) for (i,j) in zs]...,);
# (:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8, :x9, :x10, .., :x20, :x21, :x22, :x23, :x24, :x25, :x26)

struct Point{N,T}
    coords::NTuple{N,T}
end

nd(x::Point{N,T}) where {T,N} = N
Base.eltype(x::Point{N,T}) where {T,N} = T

Base.lastindex(x::Point{N,T}) where {T,N} = N
Base.lastindex(::Type{Point{N,T}}) where {T,N} = N

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
    
Base.getindex(x::Point{N,T}, idx::I) where {T,N,I<:Union{Signed,Unsigned}} = x.coords[idx]
Base.getindex(x::Point{N,T}, idxs::R) where {T,N,R<:UnitRange} = x.coords[idxs]

function Base.setindex!(pt::Point{N,T}, value::T, idx::Signed) where {T,N}
    idx == 1 && return Point(value, pt[2:end]...,)
    idx == N && return Point(pt[1:end-1]..., value)
    return Point(pt[1:(idx-1)]...,value,pt[idx+1:end]...,)
end

function Base.setindex!(pt::Point{N,T}, values::NTuple{M,T}, idxs::R) where {T,N,M,R<:UnitRange}
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

Δpt²(pt1::T) where {F, T<:Point{1,F}} = Δx1²(pt1)
Δpt²(pt1::T) where {F, T<:Point{2,F}} = Δx1²(pt1) + Δx2²(pt1) 
Δpt²(pt1::T) where {F, T<:Point{3,F}} = Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) 
Δpt²(pt1::T) where {F, T<:Point{4,F}} = Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) 
Δpt²(pt1::T) where {F, T<:Point{5,F}} = Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) + Δx5²(pt1) 
Δpt²(pt1::T) where {F, T<:Point{6,F}} = Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) + Δx5²(pt1) + Δx6²(pt1) 
Δpt²(pt1::T) where {F, T<:Point{7,F}} = Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) + Δx5²(pt1) + Δx6²(pt1) + Δx7²(pt1) 
Δpt²(pt1::T) where {F, T<:Point{8,F}} = Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) + Δx5²(pt1) + Δx6²(pt1) + Δx7²(pt1) + Δx8²(pt1)
Δpt²(pt1::T) where {F, T<:Point{9,F}} = (Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) + Δx5²(pt1) + Δx6²(pt1) + Δx7²(pt1) + Δx8²(pt1)) + Δx9²(pt1)
Δpt²(pt1::T) where {F, T<:Point{10,F}} = (Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) + Δx5²(pt1) + Δx6²(pt1) + Δx7²(pt1) + Δx8²(pt1)) + (Δx9²(pt1) + Δx10²(pt1))
Δpt²(pt1::T) where {F, T<:Point{11,F}} = (Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) + Δx5²(pt1) + Δx6²(pt1) + Δx7²(pt1) + Δx8²(pt1)) + (Δx9²(pt1) + Δx10²(pt1) + Δx11²(pt1))
Δpt²(pt1::T) where {F, T<:Point{12,F}} = (Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) + Δx5²(pt1) + Δx6²(pt1) + Δx7²(pt1) + Δx8²(pt1)) + (Δx9²(pt1) + Δx10²(pt1) + Δx11²(pt1) + Δx12²(pt1))
Δpt²(pt1::T) where {F, T<:Point{13,F}} = (Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) + Δx5²(pt1) + Δx6²(pt1) + Δx7²(pt1) + Δx8²(pt1)) + (Δx9²(pt1) + Δx10²(pt1) + Δx11²(pt1) + Δx12²(pt1) + Δx13²(pt1))
Δpt²(pt1::T) where {F, T<:Point{14,F}} = (Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) + Δx5²(pt1) + Δx6²(pt1) + Δx7²(pt1) + Δx8²(pt1)) + (Δx9²(pt1) + Δx10²(pt1) + Δx11²(pt1) + Δx12²(pt1) + Δx13²(pt1) + Δx14²(pt1))
Δpt²(pt1::T) where {F, T<:Point{15,F}} = (Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) + Δx5²(pt1) + Δx6²(pt1) + Δx7²(pt1) + Δx8²(pt1)) + (Δx9²(pt1) + Δx10²(pt1) + Δx11²(pt1) + Δx12²(pt1) + Δx13²(pt1) + Δx14²(pt1) + Δx15²(pt1))
Δpt²(pt1::T) where {F, T<:Point{16,F}} = (Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) + Δx5²(pt1) + Δx6²(pt1) + Δx7²(pt1) + Δx8²(pt1)) + (Δx9²(pt1) + Δx10²(pt1) + Δx11²(pt1) + Δx12²(pt1) + Δx13²(pt1) + Δx14²(pt1) + Δx15²(pt1) + Δx16²(pt1))
Δpt²(pt1::T) where {F, T<:Point{17,F}} = (Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) + Δx5²(pt1) + Δx6²(pt1) + Δx7²(pt1) + Δx8²(pt1)) + (Δx9²(pt1) + Δx10²(pt1) + Δx11²(pt1) + Δx12²(pt1) + Δx13²(pt1) + Δx14²(pt1) + Δx15²(pt1) + Δx16²(pt1)) + (Δx17²(pt1))
Δpt²(pt1::T) where {F, T<:Point{18,F}} = (Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) + Δx5²(pt1) + Δx6²(pt1) + Δx7²(pt1) + Δx8²(pt1)) + (Δx9²(pt1) + Δx10²(pt1) + Δx11²(pt1) + Δx12²(pt1) + Δx13²(pt1) + Δx14²(pt1) + Δx15²(pt1) + Δx16²(pt1)) + (Δx17²(pt1) + Δx18²(pt1))
Δpt²(pt1::T) where {F, T<:Point{19,F}} = (Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) + Δx5²(pt1) + Δx6²(pt1) + Δx7²(pt1) + Δx8²(pt1)) + (Δx9²(pt1) + Δx10²(pt1) + Δx11²(pt1) + Δx12²(pt1) + Δx13²(pt1) + Δx14²(pt1) + Δx15²(pt1) + Δx16²(pt1)) + (Δx17²(pt1) + Δx18²(pt1) + Δx19²(pt1))
Δpt²(pt1::T) where {F, T<:Point{20,F}} = (Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) + Δx5²(pt1) + Δx6²(pt1) + Δx7²(pt1) + Δx8²(pt1)) + (Δx9²(pt1) + Δx10²(pt1) + Δx11²(pt1) + Δx12²(pt1) + Δx13²(pt1) + Δx14²(pt1) + Δx15²(pt1) + Δx16²(pt1)) + (Δx17²(pt1) + Δx18²(pt1) + Δx19²(pt1) + Δx20²(pt1))
Δpt²(pt1::T) where {F, T<:Point{21,F}} = (Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) + Δx5²(pt1) + Δx6²(pt1) + Δx7²(pt1) + Δx8²(pt1)) + (Δx9²(pt1) + Δx10²(pt1) + Δx11²(pt1) + Δx12²(pt1) + Δx13²(pt1) + Δx14²(pt1) + Δx15²(pt1) + Δx16²(pt1)) + (Δx17²(pt1) + Δx18²(pt1) + Δx19²(pt1) + Δx20²(pt1) + Δx21²(pt1))
Δpt²(pt1::T) where {F, T<:Point{22,F}} = (Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) + Δx5²(pt1) + Δx6²(pt1) + Δx7²(pt1) + Δx8²(pt1)) + (Δx9²(pt1) + Δx10²(pt1) + Δx11²(pt1) + Δx12²(pt1) + Δx13²(pt1) + Δx14²(pt1) + Δx15²(pt1) + Δx16²(pt1)) + (Δx17²(pt1) + Δx18²(pt1) + Δx19²(pt1) + Δx20²(pt1) + Δx21²(pt1) + Δx22²(pt1))
Δpt²(pt1::T) where {F, T<:Point{23,F}} = (Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) + Δx5²(pt1) + Δx6²(pt1) + Δx7²(pt1) + Δx8²(pt1)) + (Δx9²(pt1) + Δx10²(pt1) + Δx11²(pt1) + Δx12²(pt1) + Δx13²(pt1) + Δx14²(pt1) + Δx15²(pt1) + Δx16²(pt1)) + (Δx17²(pt1) + Δx18²(pt1) + Δx19²(pt1) + Δx20²(pt1) + Δx21²(pt1) + Δx22²(pt1) + Δx23²(pt1))
Δpt²(pt1::T) where {F, T<:Point{24,F}} = (Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) + Δx5²(pt1) + Δx6²(pt1) + Δx7²(pt1) + Δx8²(pt1)) + (Δx9²(pt1) + Δx10²(pt1) + Δx11²(pt1) + Δx12²(pt1) + Δx13²(pt1) + Δx14²(pt1) + Δx15²(pt1) + Δx16²(pt1)) + (Δx17²(pt1) + Δx18²(pt1) + Δx19²(pt1) + Δx20²(pt1) + Δx21²(pt1) + Δx22²(pt1) + Δx23²(pt1) + Δx24²(pt1))
Δpt²(pt1::T) where {F, T<:Point{25,F}} = (Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) + Δx5²(pt1) + Δx6²(pt1) + Δx7²(pt1) + Δx8²(pt1)) + (Δx9²(pt1) + Δx10²(pt1) + Δx11²(pt1) + Δx12²(pt1) + Δx13²(pt1) + Δx14²(pt1) + Δx15²(pt1) + Δx16²(pt1)) + (Δx17²(pt1) + Δx18²(pt1) + Δx19²(pt1) + Δx20²(pt1) + Δx21²(pt1) + Δx22²(pt1) + Δx23²(pt1) + Δx24²(pt1)) + Δx25²(pt1)
Δpt²(pt1::T) where {F, T<:Point{26,F}} = (Δx1²(pt1) + Δx2²(pt1) + Δx3²(pt1) + Δx4²(pt1) + Δx5²(pt1) + Δx6²(pt1) + Δx7²(pt1) + Δx8²(pt1)) + (Δx9²(pt1) + Δx10²(pt1) + Δx11²(pt1) + Δx12²(pt1) + Δx13²(pt1) + Δx14²(pt1) + Δx15²(pt1) + Δx16²(pt1)) + (Δx17²(pt1) + Δx18²(pt1) + Δx19²(pt1) + Δx20²(pt1) + Δx21²(pt1) + Δx22²(pt1) + Δx23²(pt1) + Δx24²(pt1)) + (Δx25²(pt1) + Δx26²(pt1))

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

Δx1(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x1(nt2) - x1(nt1)
Δx2(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x2(nt2) - x2(nt1)
Δx3(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x3(nt2) - x3(nt1)
Δx4(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x4(nt2) - x4(nt1)
Δx5(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x5(nt2) - x5(nt1)
Δx6(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x6(nt2) - x6(nt1)
Δx7(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x7(nt2) - x7(nt1)
Δx8(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x8(nt2) - x8(nt1)
Δx9(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x9(nt2) - x9(nt1)
Δx10(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x10(nt2) - x10(nt1)
Δx11(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x11(nt2) - x11(nt1)
Δx12(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x12(nt2) - x12(nt1)
Δx13(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x13(nt2) - x13(nt1)
Δx14(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x14(nt2) - x14(nt1)
Δx15(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x15(nt2) - x15(nt1)
Δx16(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x16(nt2) - x16(nt1)
Δx17(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x17(nt2) - x17(nt1)
Δx18(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x18(nt2) - x18(nt1)
Δx19(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x19(nt2) - x19(nt1)
Δx20(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x20(nt2) - x20(nt1)
Δx21(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x21(nt2) - x21(nt1)
Δx22(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x22(nt2) - x22(nt1)
Δx23(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x23(nt2) - x23(nt1)
Δx24(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x24(nt2) - x24(nt1)
Δx25(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x25(nt2) - x25(nt1)
Δx26(nt1::T, nt2::T) where {F, T<:Point{N,F}} = x26(nt2) - x26(nt1)

# squared separation for coordinate axes, unoriented

Δx1²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx1(nt1::T, nt2::T); d*d end
Δx2²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx2(nt1::T, nt2::T); d*d end
Δx3²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx3(nt1::T, nt2::T); d*d end
Δx4²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx4(nt1::T, nt2::T); d*d end
Δx5²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx5(nt1::T, nt2::T); d*d end
Δx6²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx6(nt1::T, nt2::T); d*d end
Δx7²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx7(nt1::T, nt2::T); d*d end
Δx8²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx8(nt1::T, nt2::T); d*d end
Δx9²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx9(nt1::T, nt2::T); d*d end
Δx10²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx10(nt1::T, nt2::T); d*d end
Δx11²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx11(nt1::T, nt2::T); d*d end
Δx12²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx12(nt1::T, nt2::T); d*d end
Δx13²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx13(nt1::T, nt2::T); d*d end
Δx14²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx14(nt1::T, nt2::T); d*d end
Δx15²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx15(nt1::T, nt2::T); d*d end
Δx16²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx16(nt1::T, nt2::T); d*d end
Δx17²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx17(nt1::T, nt2::T); d*d end
Δx18²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx18(nt1::T, nt2::T); d*d end
Δx19²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx19(nt1::T, nt2::T); d*d end
Δx20²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx20(nt1::T, nt2::T); d*d end
Δx21²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx21(nt1::T, nt2::T); d*d end
Δx22²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx22(nt1::T, nt2::T); d*d end
Δx23²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx23(nt1::T, nt2::T); d*d end
Δx24²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx24(nt1::T, nt2::T); d*d end
Δx25²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx25(nt1::T, nt2::T); d*d end
Δx26²(nt1::T, nt2::T) where {F, T<:Point{N,F}} = let d = Δx26(nt1::T, nt2::T); d*d end

# squared interpoint distance (norm2)

Δpt²(pt1::T, pt2::T) where {F, T<:Point{1,F}} = Δx1²(pt1, pt2)
Δpt²(pt1::T, pt2::T) where {F, T<:Point{2,F}} = Δx1²(pt1, pt2) + Δx2²(pt1, pt2) 
Δpt²(pt1::T, pt2::T) where {F, T<:Point{3,F}} = Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) 
Δpt²(pt1::T, pt2::T) where {F, T<:Point{4,F}} = Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) 
Δpt²(pt1::T, pt2::T) where {F, T<:Point{5,F}} = Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) + Δx5²(pt1, pt2) 
Δpt²(pt1::T, pt2::T) where {F, T<:Point{6,F}} = Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) + Δx5²(pt1, pt2) + Δx6²(pt1, pt2) 
Δpt²(pt1::T, pt2::T) where {F, T<:Point{7,F}} = Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) + Δx5²(pt1, pt2) + Δx6²(pt1, pt2) + Δx7²(pt1, pt2) 
Δpt²(pt1::T, pt2::T) where {F, T<:Point{8,F}} = Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) + Δx5²(pt1, pt2) + Δx6²(pt1, pt2) + Δx7²(pt1, pt2) + Δx8²(pt1, pt2)
Δpt²(pt1::T, pt2::T) where {F, T<:Point{9,F}} = (Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) + Δx5²(pt1, pt2) + Δx6²(pt1, pt2) + Δx7²(pt1, pt2) + Δx8²(pt1, pt2)) + Δx9²(pt1, pt2)
Δpt²(pt1::T, pt2::T) where {F, T<:Point{10,F}} = (Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) + Δx5²(pt1, pt2) + Δx6²(pt1, pt2) + Δx7²(pt1, pt2) + Δx8²(pt1, pt2)) + (Δx9²(pt1, pt2) + Δx10²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{11,F}} = (Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) + Δx5²(pt1, pt2) + Δx6²(pt1, pt2) + Δx7²(pt1, pt2) + Δx8²(pt1, pt2)) + (Δx9²(pt1, pt2) + Δx10²(pt1, pt2) + Δx11²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{12,F}} = (Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) + Δx5²(pt1, pt2) + Δx6²(pt1, pt2) + Δx7²(pt1, pt2) + Δx8²(pt1, pt2)) + (Δx9²(pt1, pt2) + Δx10²(pt1, pt2) + Δx11²(pt1, pt2) + Δx12²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{13,F}} = (Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) + Δx5²(pt1, pt2) + Δx6²(pt1, pt2) + Δx7²(pt1, pt2) + Δx8²(pt1, pt2)) + (Δx9²(pt1, pt2) + Δx10²(pt1, pt2) + Δx11²(pt1, pt2) + Δx12²(pt1, pt2) + Δx13²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{14,F}} = (Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) + Δx5²(pt1, pt2) + Δx6²(pt1, pt2) + Δx7²(pt1, pt2) + Δx8²(pt1, pt2)) + (Δx9²(pt1, pt2) + Δx10²(pt1, pt2) + Δx11²(pt1, pt2) + Δx12²(pt1, pt2) + Δx13²(pt1, pt2) + Δx14²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{15,F}} = (Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) + Δx5²(pt1, pt2) + Δx6²(pt1, pt2) + Δx7²(pt1, pt2) + Δx8²(pt1, pt2)) + (Δx9²(pt1, pt2) + Δx10²(pt1, pt2) + Δx11²(pt1, pt2) + Δx12²(pt1, pt2) + Δx13²(pt1, pt2) + Δx14²(pt1, pt2) + Δx15²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{16,F}} = (Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) + Δx5²(pt1, pt2) + Δx6²(pt1, pt2) + Δx7²(pt1, pt2) + Δx8²(pt1, pt2)) + (Δx9²(pt1, pt2) + Δx10²(pt1, pt2) + Δx11²(pt1, pt2) + Δx12²(pt1, pt2) + Δx13²(pt1, pt2) + Δx14²(pt1, pt2) + Δx15²(pt1, pt2) + Δx16²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{17,F}} = (Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) + Δx5²(pt1, pt2) + Δx6²(pt1, pt2) + Δx7²(pt1, pt2) + Δx8²(pt1, pt2)) + (Δx9²(pt1, pt2) + Δx10²(pt1, pt2) + Δx11²(pt1, pt2) + Δx12²(pt1, pt2) + Δx13²(pt1, pt2) + Δx14²(pt1, pt2) + Δx15²(pt1, pt2) + Δx16²(pt1, pt2)) + (Δx17²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{18,F}} = (Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) + Δx5²(pt1, pt2) + Δx6²(pt1, pt2) + Δx7²(pt1, pt2) + Δx8²(pt1, pt2)) + (Δx9²(pt1, pt2) + Δx10²(pt1, pt2) + Δx11²(pt1, pt2) + Δx12²(pt1, pt2) + Δx13²(pt1, pt2) + Δx14²(pt1, pt2) + Δx15²(pt1, pt2) + Δx16²(pt1, pt2)) + (Δx17²(pt1, pt2) + Δx18²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{19,F}} = (Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) + Δx5²(pt1, pt2) + Δx6²(pt1, pt2) + Δx7²(pt1, pt2) + Δx8²(pt1, pt2)) + (Δx9²(pt1, pt2) + Δx10²(pt1, pt2) + Δx11²(pt1, pt2) + Δx12²(pt1, pt2) + Δx13²(pt1, pt2) + Δx14²(pt1, pt2) + Δx15²(pt1, pt2) + Δx16²(pt1, pt2)) + (Δx17²(pt1, pt2) + Δx18²(pt1, pt2) + Δx19²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{20,F}} = (Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) + Δx5²(pt1, pt2) + Δx6²(pt1, pt2) + Δx7²(pt1, pt2) + Δx8²(pt1, pt2)) + (Δx9²(pt1, pt2) + Δx10²(pt1, pt2) + Δx11²(pt1, pt2) + Δx12²(pt1, pt2) + Δx13²(pt1, pt2) + Δx14²(pt1, pt2) + Δx15²(pt1, pt2) + Δx16²(pt1, pt2)) + (Δx17²(pt1, pt2) + Δx18²(pt1, pt2) + Δx19²(pt1, pt2) + Δx20²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{21,F}} = (Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) + Δx5²(pt1, pt2) + Δx6²(pt1, pt2) + Δx7²(pt1, pt2) + Δx8²(pt1, pt2)) + (Δx9²(pt1, pt2) + Δx10²(pt1, pt2) + Δx11²(pt1, pt2) + Δx12²(pt1, pt2) + Δx13²(pt1, pt2) + Δx14²(pt1, pt2) + Δx15²(pt1, pt2) + Δx16²(pt1, pt2)) + (Δx17²(pt1, pt2) + Δx18²(pt1, pt2) + Δx19²(pt1, pt2) + Δx20²(pt1, pt2) + Δx21²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{22,F}} = (Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) + Δx5²(pt1, pt2) + Δx6²(pt1, pt2) + Δx7²(pt1, pt2) + Δx8²(pt1, pt2)) + (Δx9²(pt1, pt2) + Δx10²(pt1, pt2) + Δx11²(pt1, pt2) + Δx12²(pt1, pt2) + Δx13²(pt1, pt2) + Δx14²(pt1, pt2) + Δx15²(pt1, pt2) + Δx16²(pt1, pt2)) + (Δx17²(pt1, pt2) + Δx18²(pt1, pt2) + Δx19²(pt1, pt2) + Δx20²(pt1, pt2) + Δx21²(pt1, pt2) + Δx22²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{23,F}} = (Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) + Δx5²(pt1, pt2) + Δx6²(pt1, pt2) + Δx7²(pt1, pt2) + Δx8²(pt1, pt2)) + (Δx9²(pt1, pt2) + Δx10²(pt1, pt2) + Δx11²(pt1, pt2) + Δx12²(pt1, pt2) + Δx13²(pt1, pt2) + Δx14²(pt1, pt2) + Δx15²(pt1, pt2) + Δx16²(pt1, pt2)) + (Δx17²(pt1, pt2) + Δx18²(pt1, pt2) + Δx19²(pt1, pt2) + Δx20²(pt1, pt2) + Δx21²(pt1, pt2) + Δx22²(pt1, pt2) + Δx23²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{24,F}} = (Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) + Δx5²(pt1, pt2) + Δx6²(pt1, pt2) + Δx7²(pt1, pt2) + Δx8²(pt1, pt2)) + (Δx9²(pt1, pt2) + Δx10²(pt1, pt2) + Δx11²(pt1, pt2) + Δx12²(pt1, pt2) + Δx13²(pt1, pt2) + Δx14²(pt1, pt2) + Δx15²(pt1, pt2) + Δx16²(pt1, pt2)) + (Δx17²(pt1, pt2) + Δx18²(pt1, pt2) + Δx19²(pt1, pt2) + Δx20²(pt1, pt2) + Δx21²(pt1, pt2) + Δx22²(pt1, pt2) + Δx23²(pt1, pt2) + Δx24²(pt1, pt2))
Δpt²(pt1::T, pt2::T) where {F, T<:Point{25,F}} = (Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) + Δx5²(pt1, pt2) + Δx6²(pt1, pt2) + Δx7²(pt1, pt2) + Δx8²(pt1, pt2)) + (Δx9²(pt1, pt2) + Δx10²(pt1, pt2) + Δx11²(pt1, pt2) + Δx12²(pt1, pt2) + Δx13²(pt1, pt2) + Δx14²(pt1, pt2) + Δx15²(pt1, pt2) + Δx16²(pt1, pt2)) + (Δx17²(pt1, pt2) + Δx18²(pt1, pt2) + Δx19²(pt1, pt2) + Δx20²(pt1, pt2) + Δx21²(pt1, pt2) + Δx22²(pt1, pt2) + Δx23²(pt1, pt2) + Δx24²(pt1, pt2)) + Δx25²(pt1, pt2)
Δpt²(pt1::T, pt2::T) where {F, T<:Point{26,F}} = (Δx1²(pt1, pt2) + Δx2²(pt1, pt2) + Δx3²(pt1, pt2) + Δx4²(pt1, pt2) + Δx5²(pt1, pt2) + Δx6²(pt1, pt2) + Δx7²(pt1, pt2) + Δx8²(pt1, pt2)) + (Δx9²(pt1, pt2) + Δx10²(pt1, pt2) + Δx11²(pt1, pt2) + Δx12²(pt1, pt2) + Δx13²(pt1, pt2) + Δx14²(pt1, pt2) + Δx15²(pt1, pt2) + Δx16²(pt1, pt2)) + (Δx17²(pt1, pt2) + Δx18²(pt1, pt2) + Δx19²(pt1, pt2) + Δx20²(pt1, pt2) + Δx21²(pt1, pt2) + Δx22²(pt1, pt2) + Δx23²(pt1, pt2) + Δx24²(pt1, pt2)) + (Δx25²(pt1, pt2) + Δx26²(pt1, pt2))


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

