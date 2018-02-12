if nworkers() > 1
    rmprocs(workers())
end

if Base.JLOptions().code_coverage == 1
    addprocs(Sys.CPU_CORES, exeflags = ["--code-coverage=user", "--inline=no", "--check-bounds=yes"])
else
    addprocs(Sys.CPU_CORES, exeflags = "--check-bounds=yes")
end
@show nprocs()

using InvariantDistribution
using SimplexSplitting
using SimplexIntersection
using Distributions
using Base.Test

include("test_mm_sparse.jl")
#include("test_markovmatrix.jl")
