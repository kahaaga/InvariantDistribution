module InvariantDistribution

using SimplexSplitting
using SimplexIntersection

include("markovmatrix.jl")
include("statespace.jl")
include("invariantdist.jl")
include("jointdist.jl")
include("rowindexin.jl")
include("transferentropy.jl")
include("te_test.jl")

export markovmatrix, markovmatrix_parallel, invariantdist, StateSpace, get_nonempty_bins,
    jointdist, marginaldists, indexin_rows, transferentropy, te_test

# package code goes here

end # module
