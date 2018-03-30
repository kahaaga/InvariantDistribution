__precompile__(true)

module InvariantDistribution

using Simplices, SimplexSplitting

include("invariant_maps.jl")
include("markovmatrix.jl")
include("statespace.jl")
include("example_embeddings.jl")
include("jointdist.jl")
include("invdist.jl")
include("invariantize_embedding.jl")
include("rowindexin.jl")
include("transferentropy.jl")
include("te_test.jl")
include("mm_sparse.jl")
include("mm_sparse_parallel.jl")
include("mm_parallel.jl")
include("mm_discrete_dense.jl")

export markovmatrix,
    markovmatrixp,
    estimate_invdist,
    StateSpace,
    get_nonempty_bins,
    jointdist,
    marginaldists,
    indexin_rows,
    invariant_gaussian_embedding,
    transferentropy,
    te_test,
    mm_sparse,
    mm_sparse_parallel,
    mm_p,
    mm_discrete_dense,
    invariantize_embedding,
    invariant_bullmap_on_cube

# Run some examples to trigger precompilation
t = SimplexSplitting.triang_from_embedding(SimplexSplitting.Embedding(InvariantDistribution.invariant_gaussian_embedding(npts = 10, covariance = 0.3, tau = 1)))
mm_discrete_dense(t, n_randpts = 10)
mm_p(t)

end # module
