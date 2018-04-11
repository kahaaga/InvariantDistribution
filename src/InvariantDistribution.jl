__precompile__(true)

module InvariantDistribution

include("invariant_maps.jl")
include("statespace.jl")
include("example_embeddings.jl")
include("left_eigenvector.jl")
include("invariantize_embedding.jl")
include("mm_exact/mm_sparse.jl")
include("mm_exact/mm_sparse_parallel.jl")
include("mm_exact/mm_parallel.jl")
include("mm_approx/mm_dd.jl")
include("mm_approx/mm_dd3.jl")
include("mm_approx/mm_dd2.jl")
include("mm_approx/mm_dd4.jl")

export markovmatrix,
    markovmatrixp,
    left_eigenvector,
    StateSpace,
    indexin_rows,
    invariant_gaussian_embedding,
    te_test,
    mm_sparse,
    mm_sparse_parallel,
    mm_p,
    invariantize_embedding,
    invariant_bullmap_on_cube,
    prepare_mm_dd,
	mm_dd,
    mm_dd3,
    mm_dd2,
    mm_dd3,
    mm_dd4

t = triang_from_embedding(Embedding(invariant_gaussian_embedding(npts = 14, cov = 0.5, tau = 1)))
mm_dd4(t)

end # module
