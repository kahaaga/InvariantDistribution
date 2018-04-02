# if nworkers() > 1
#     rmprocs(workers())
# end
#
# if Base.JLOptions().code_coverage == 1
#     addprocs(Sys.CPU_CORES - 3 , exeflags = ["--code-coverage=user", "--inline=no", "--check-bounds=yes"])
# else
#     addprocs(Sys.CPU_CORES - 3, exeflags = "--check-bounds=yes")
# end
# @show nprocs()

using Simplices
using SimplexSplitting
using InvariantDistribution
using Distributions
using Base.Test

include("mm_dd_all.jl")
include("invdist.jl")
