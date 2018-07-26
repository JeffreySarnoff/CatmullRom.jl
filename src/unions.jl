const TupOfValues{L,F} = NTuple{L,F} where {L,F}
const VecOfValues{F}   = Vector{F} where {F}
const ValueSeq{L,F}    = Union{TupOfValues{L,F}, VecOfValues{F}} where {L,F}

const TupAsPoint{D,R} = NTuple{D,R} where {D,R}
const VecAsPoint{R}   = Vector{R}   where {R}

const TupOfTupPoints{M,D,R} = NTuple{M, TupAsPoint{D,R}} where {M,D,R}
const TupOfVecPoints{M,R}   = NTuple{M, VecAsPoint{R}}   where {M,R}
const VecOfTupPoints{D,R}   = Vector{TupAsPoint{D,R}}    where {D,R}
const VecOfVecPoints{R}     = Vector{VecAsPoint{R}}      where {R}

const PointSeq{M,D,R}  = Union{TupOfTupPoints{M,D,R}, TupOfVecPoints{M,R}, 
                               VecOfTupPoints{D,R}, VecOfVecPoints{R}} where {M,D,R}
