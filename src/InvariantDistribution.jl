__precompile__()

module InvariantDistribution

using Simplices
using SimplexSplitting

include("invariant_maps.jl")
include("markovmatrix.jl")
include("statespace.jl")
include("invariantdist.jl")
include("jointdist.jl")
include("rowindexin.jl")
include("transferentropy.jl")
include("draw_invariant_embedding.jl")
include("te_test.jl")
include("mm_sparse.jl")
include("mm_sparse_parallel.jl")
include("mm_parallel.jl")
include("mm_discrete_dense.jl")

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
    mm_discrete_dense_test,
    mm_discrete_dense,
    invariant_gaussian_embedding,
    invariantize_embedding,
    invariant_bullmap_on_cube,
    test_mm_discrete_map

# package code goes here

end # module
