const Open   = :Open
const Closed = :Closed

const Linear    = :Linear
const Quadratic = :Quadratic
const Thiele3   = :Thiele3
const Omit      = :Omit

# stubs
linear()    = ()
quadratic() = ()
thiele3()   = ()

const endpointfn = Dict([:Linear => linear, :Quadratic => quadratic, :Thiele3 => thiele3])
