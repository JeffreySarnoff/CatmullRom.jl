"""
    ProtoNT(names::NTuple{N,Symbol})
    
Generate a NamedTuple Prototype that labels precipient values using `names`.
"""
macro ProtoNT(names)
   :(NamedTuple{$names})
end


"""
    Point in 1D, 2D, 3D, 4D

   a2Dpoint = PT2D( (1.0, 2.0) )
   a2Dpoint = @PT2D( 1.0, 2.0 )
"""

PT1D = @ProtoNT( (:x,) )
PT2D = @ProtoNT( (:x, :y) )
PT3D = @ProtoNT( (:x, :y, :z) )
PT4D = @ProtoNT( (:x, :y, :z, :t) )

macro PT1D(x)
    :(PT1D(($x,)))
end

macro PT2D(x, y)
    :(PT2D(($x, $y)))
end

macro PT3D(x, y, z)
    :(PT3D(($x, $y, $z)))
end

macro PT4D(x, y, z, t)
    :(PT4D(($x, $y, $z, $t)))
end


# retrieve the coordinates of a Point as a Tuple

coords(nt::PT1D) = values(nt)
coords(nt::PT2D) = values(nt)
coords(nt::PT3D) = values(nt)
coords(nt::PT4D) = values(nt)

# retrieve a coordinate from a Point as a Number

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


# separation for coordinate axes, oriented

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
