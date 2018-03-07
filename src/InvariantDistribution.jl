__precompile__()

module InvariantDistribution

using Simplices
using SimplexSplitting

include("markovmatrix.jl")
include("statespace.jl")
include("invariantdist.jl")
include("jointdist.jl")
include("rowindexin.jl")
include("transferentropy.jl")
include("te_test.jl")
include("mm_sparse.jl")
include("mm_sparse_parallel.jl")
include("mm_parallel.jl")
include("draw_invariant_embedding.jl")

export markovmatrix,
    markovmatrixp,
    invariantdist,
    StateSpace,
    get_nonempty_bins,
    jointdist,
    marginaldists,
    indexin_rows,
    transferentropy,
    te_test,
    mm_sparse,
    mm_sparse_parallel,
    mm_p,
    invariant_gaussian_embedding,
    invariantize_embedding

# package code goes here

end # module
