VERSION < v"0.6" && warn("The package has not been tested on Julia versions < 0.6")
using SimplexSplitting
using BenchmarkTools
using InvariantDistribution

# Define a parent BenchmarkGroup to contain our suite
const suite = BenchmarkGroup()

#########
# Setup
#########

# Create a triangulation to work with
t = SimplexSplitting.triang_from_embedding(SimplexSplitting.Embedding(InvariantDistribution.invariant_gaussian_embedding(npts = 30, covariance = 0.3, tau = 1)))

# Add some child groups to our benchmark suite.
suite["markovmatrix"] = BenchmarkGroup(["discrete", "exact"])

#########
# Benchmark
#########
#suite["markovmatrix"]["exact"] = @benchmarkable mm_p(t)
suite["markovmatrix"]["discrete"] = @benchmarkable mm_discrete_dense(t, n_randpts = 100)

# If a cache of tuned parameters already exists, use it, otherwise, tune and cache
# the benchmark parameters. Reusing cached parameters is faster and more reliable
# than re-tuning `suite` every time the file is included.
paramspath = joinpath(dirname(@__FILE__), "params.json")

if isfile(paramspath)
    loadparams!(suite, BenchmarkTools.load(paramspath)[1], :evals);
else
    tune!(suite)
    BenchmarkTools.save(paramspath, params(suite));
end
