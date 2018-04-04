installed = Pkg.installed()
if !("Simplices" in keys(installed))
    Pkg.clone("https://github.com/kahaaga/Simplices.jl")
else
    using Simplices
end

using Simplices
using SimplexSplitting
using InvariantDistribution
using Distributions
using Base.Test

include("mm_dd_all.jl")
include("invdist.jl")
