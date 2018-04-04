__precompile__(true)

module InvariantDistribution

using Simplices, SimplexSplitting, Parameters,
    StaticArrays, Distributions, InplaceOps

installed = Pkg.installed()
if !("Pycall" in keys(installed))
    Pkg.add("PyCall")
    ENV["PYTHON"]= ""; Pkg.build("PyCall")
end

if !("Conda" in keys(installed))
    Pkg.add("Conda")
    using Conda; Conda.add("scipy")
end

if !("Simplices" in keys(installed))
    Pkg.clone("https://github.com/kahaaga/Simplices.jl")
else
    using Simplices
end


include("invariant_maps.jl")
include("statespace.jl")
include("example_embeddings.jl")
include("invdist.jl")
include("invariantize_embedding.jl")
include("mm_exact/mm_sparse.jl")
include("mm_exact/mm_sparse_parallel.jl")
include("mm_exact/mm_parallel.jl")
include("mm_approx/mm_dd.jl")
include("mm_approx/mm_dd2.jl")
include("mm_approx/mm_dd4.jl")

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
    prepare_mm_dd,
	mm_dd,
    mm_dd2,
    mm_dd3,
    mm_dd4

t = triang_from_embedding(Embedding(invariant_gaussian_embedding(npts = 8, covariance = 0.5, tau = 1)))
mm_dd4(t)

end # module
