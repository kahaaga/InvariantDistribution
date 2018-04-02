__precompile__(true)

module InvariantDistribution

using Simplices, SimplexSplitting

include("invariant_maps.jl")
include("statespace.jl")
include("example_embeddings.jl")
include("invdist.jl")
include("invariantize_embedding.jl")
include("mm_sparse.jl")
include("mm_sparse_parallel.jl")
include("mm_parallel.jl")
include("mm_dd.jl")
include("mm_dd2.jl")

export markovmatrix,
    markovmatrixp,
    estimate_invdist,
    StateSpace,
    indexin_rows,
    invariant_gaussian_embedding,
    te_test,
    mm_sparse,
    mm_sparse_parallel,
    mm_p,
    invariantize_embedding,
    invariant_bullmap_on_cube,
	mm_dd, mm_dd2, mm_dd3, get_coeffs, prepare_mm_dd


function m1(t::Triangulation; n_pts = 100, sample_randomly = false)
	mm_dd(t, n_randpts = n_pts, sample_randomly = sample_randomly)
end

function m2(t::Triangulation; n_pts = 100, sample_randomly = false)
	mm_dd2(t, n_randpts = n_pts, sample_randomly = sample_randomly)
end

function m3(t::Triangulation; n_pts = 100, sample_randomly = false)
	dim = size(t.points, 2)
	convex_coeffs = get_coeffs(dim, n_pts, sample_randomly)
	simplices, imsimplices = prepare_mm_dd(t)
	mm_dd3(t, simplices, convex_coeffs)
end

# Run some examples to trigger precompilation
t = triang_from_embedding(Embedding(invariant_gaussian_embedding(npts = 8, covariance = 0.5, tau = 1)))
m1(t, n_pts = 10)
m2(t, n_pts = 10)
m3(t, n_pts = 10)

end # module
