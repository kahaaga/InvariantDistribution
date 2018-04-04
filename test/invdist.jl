using SimplexSplitting, InvariantDistribution

@testset "InvDist" begin
    e_ex = Embedding(invariant_gaussian_embedding(npts = 25))
    t_ex = triang_from_embedding(e_ex)
    mm = mm_dd(t_ex)
    @test all(sum(mm, 2) .≈ 1)

    invdist = estimate_invdist(mm)
    @test sum(invdist.dist) ≈ 1
end
