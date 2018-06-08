"""
     Prototypes for NamedTuples
     
     Specify the names without specifying their value's types.
     
     The names are given as a tuple of symbols.
     
     ```julia
     using Dates
     
     field_names  = (:city, :timezone, :utcplus)
     TimeZoneInfo = ProtoNT( field_names )
     
     city     = "Boston"
     timezone = "America/NewYork"
     utcplus  = (Hour(-5), Hour(-4))
     
     boston_zoneinfo = TimeZoneInfo( (city, timezone, utcplus) ) 
     ```
     
""" ProtoNT

ProtoNT(z::Symbol) = NamedTuple{(z,)}
ProtoNT(z::Tuple)  = NamedTuple{(map(Symbol, z)...,)}
ProtoNT(z::Vector) = NamedTuple{(map(Symbol, z)...,)}
ProtoNT(z::Symbol, zs...) = NamedTuple{(z, map(Symbol, zs)...,)}
ProtoNT(z::AbstractString) = NamedTuple{(Symbol(z),)}
ProtoNT(z::AbstractString, zs...) = NamedTuple{(Symbol(z), map(Symbol, zs)...,)}

"""
    Point in 1D, 2D, 3D, 4D
"""

PT1D = ProtoNT( :x )
PT2D = ProtoNT( :x, :y )
PT3D = ProtoNT( :x, :y, :z )
PT4D = ProtoNT( :x, :y, :z, :t )

xcoord(nt::PT1D) = nt.x
xcoord(nt::PT2D) = nt.x
ycoord(nt::PT2D) = nt.y
xcoord(nt::PT3D) = nt.x
ycoord(nt::PT3D) = nt.y
zcoord(nt::PT3D) = nt.z
xcoord(nt::PT4D) = nt.x
ycoord(nt::PT4D) = nt.y
zcoord(nt::PT4D) = nt.z
tcoord(nt::PT4D) = nt.t

# distance, separation
dxcoord(nt1::T, nt2::T) where {T<:NamedTuple} = xcoord(nt2) - xcoord(nt1)
dycoord(nt1::T, nt2::T) where {T<:NamedTuple} = ycoord(nt2) - ycoord(nt1)
dzcoord(nt1::T, nt2::T) where {T<:NamedTuple} = zcoord(nt2) - zcoord(nt1)
dtcoord(nt1::T, nt2::T) where {T<:NamedTuple} = tcoord(nt2) - tcoord(nt1)

# distance squared (norm2)
dxcoord2(nt1::T, nt2::T) where {T<:NamedTuple} = let d = dxcoord(nt1, nt2); d*d; end
dycoord2(nt1::T, nt2::T) where {T<:NamedTuple} = let d = dycoord(nt1, nt2); d*d; end
dzcoord2(nt1::T, nt2::T) where {T<:NamedTuple} = let d = dzcoord(nt1, nt2); d*d; end
dtcoord2(nt1::T, nt2::T) where {T<:NamedTuple} = let d = dtcoord(nt1, nt2); d*d; end

