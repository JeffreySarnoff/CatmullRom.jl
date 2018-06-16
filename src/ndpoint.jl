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


#=
# squared separation for coordinate axes, unoriented 

Δxcoord2(nt1::T, nt2::T) where {T<:NamedTuple} = let d = Δxcoord(nt1, nt2); d*d; end
Δycoord2(nt1::T, nt2::T) where {T<:NamedTuple} = let d = Δycoord(nt1, nt2); d*d; end
Δzcoord2(nt1::T, nt2::T) where {T<:NamedTuple} = let d = Δzcoord(nt1, nt2); d*d; end
Δtcoord2(nt1::T, nt2::T) where {T<:NamedTuple} = let d = Δtcoord(nt1, nt2); d*d; end

# squared interpoint distance (norm2)

Δpoint2(pt1::T, pt2::T) where {T<:PT1D} = Δxcoord2(pt1, pt2)
Δpoint2(pt1::T, pt2::T) where {T<:PT2D} = Δxcoord2(pt1, pt2) + Δycoord2(pt1, pt2)
Δpoint2(pt1::T, pt2::T) where {T<:PT3D} = Δxcoord2(pt1, pt2) + Δycoord2(pt1, pt2) + Δzcoord(pt1, pt2)
Δpoint2(pt1::T, pt2::T) where {T<:PT4D} = Δxcoord2(pt1, pt2) + Δycoord2(pt1, pt2) + Δzcoord(pt1, pt2) + Δtcoord(pt1, pt2)

const dpoint2 = Δpoint2

# interpoint distance

Δpoint(nt1::T, nt2::T) where {T<:NamedTuple} = sqrt(Δpoint2(nt1, nt2))
Δpoint(nt1::T, nt2::T) where {T<:NamedTuple} = sqrt(Δpoint2(nt1, nt2))
Δpoint(nt1::T, nt2::T) where {T<:NamedTuple} = sqrt(Δpoint2(nt1, nt2))
Δpoint(nt1::T, nt2::T) where {T<:NamedTuple} = sqrt(Δpoint2(nt1, nt2))

const dpoint = Δpoint
=#
