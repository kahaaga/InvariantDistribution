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
include("mm_sparse.jl")
include("mm_parallel.jl")
include("draw_invariant_embedding.jl")

export markovmatrix, markovmatrixp, invariantdist, StateSpace, get_nonempty_bins,
    jointdist, marginaldists, indexin_rows, transferentropy, te_test, mm_sparse,
    invariant_embedding, mm_p

# package code goes here

end # module
