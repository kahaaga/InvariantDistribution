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
    mm_discrete_dense,
    invariant_gaussian_embedding,
    invariantize_embedding,
    invariant_bullmap_on_cube
# Run some examples to trigger precompilation
t = SimplexSplitting.triang_from_embedding(SimplexSplitting.Embedding(InvariantDistribution.invariant_gaussian_embedding(npts = 10, covariance = 0.3, tau = 1)))
mm_discrete_dense(t, n_randpts = 10)
mm_p(t)

end # module
