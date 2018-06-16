# form points in 1D..26D coordinate space

# zs = zip(split(repeat("x",26),""), string.(collect(1:26)));
# coord_symbols = ([Symbol(string(i,j)) for (i,j) in zs]...,);
# (:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8, :x9, :x10, .., :x20, :x21, :x22, :x23, :x24, :x25, :x26)

struct Point{N,T}
    coords::NTuple{N,T}
end

nd(x::Point{N,T}) where {T,N} = N
Base.eltype(x::Point{N,T}) where {T,N} = T

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
