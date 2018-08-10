const Linear    = :Linear
const Quadratic = :Quadratic
const Thiele3   = :Thiele3
const Thiele4   = :Thiele4
const Omit      = :Omit

# stubs
linear()    = ()
quadratic() = ()
thiele3()   = ()
thiele4()   = ()

const endpointfn = Dict([:Linear => linear, :Quadratic => quadratic, :Thiele3 => thiele3, :Thiele4 => thiele4])
