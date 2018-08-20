#=
A Point is represented as a Tuple of coordinate values.
All coordinates of a Point are of some shared type.

Points are represented as a Point Vector.
Points are an ordered sequence of Point occurances.
The order determines the interpolatory curve traced in space.
The sequence which may overlap at the very start/end to form a closed curve.

Interpolants are 0..1 (inclusive) ordered sequence,
    represented as a tuple of floats.
=#

const OnePoint{N,T} = NTuple{N,T} where {N, T<:Number}
const Points{N,T} = Vector{OnePoint{N,T}} where {N, T<:Number}

const Seq01{N,T} = NTuple{N,T} where {N, T<:AbstractFloat}

# extrapolation methods to include endpoints

const Linear    = :Linear
const Quadratic = :Quadratic
const Thiele3   = :Thiele3
const Thiele4   = :Thiele4
const Omit      = :Omit
