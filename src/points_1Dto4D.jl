#=
    Tuple of Values from NamedTuples

tuple_of_values(namedtuple::T) where T<:NamedTuple = Base.values(namedtuple)
=#

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

# oriented separation for coordinate axes
Δxcoord(nt1::T, nt2::T) where {T<:NamedTuple} = xcoord(nt2) - xcoord(nt1)
Δycoord(nt1::T, nt2::T) where {T<:NamedTuple} = ycoord(nt2) - ycoord(nt1)
Δzcoord(nt1::T, nt2::T) where {T<:NamedTuple} = zcoord(nt2) - zcoord(nt1)
Δtcoord(nt1::T, nt2::T) where {T<:NamedTuple} = tcoord(nt2) - tcoord(nt1)

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

const dpoint  = Δpoint
