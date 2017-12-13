module InvariantDistribution

using SimplexSplitting
using SimplexIntersection

include("markovmatrix.jl")
include("statespace.jl")
include("invariantdist.jl")
include("jointdist.jl")
include("rowindexin.jl")

export markovmatrix, markovmatrix_parallel, invariantdist, StateSpace, get_nonempty_bins,
    jointdist, marginaldists, indexin_rows

# package code goes here

end # module
