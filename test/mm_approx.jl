using SimplexSplitting, InvariantDistribution

function m1(t::Triangulation; n_randpts = 100, sample_randomly = false)
	mm_dd(t, n_randpts = n_randpts, sample_randomly = sample_randomly)
end

function m2(t::Triangulation; n_randpts = 100, sample_randomly = false)
	mm_dd2(t, n_randpts = n_randpts, sample_randomly = sample_randomly)
end

function m3(t::Triangulation; n_randpts = 100, sample_randomly = false)
	dim = size(t.points, 2)
	convex_coeffs = get_coeffs(dim, n_randpts, sample_randomly)
	simplices, imsimplices = prepare_mm_dd(t)
	mm_dd3(t, simplices, convex_coeffs)
end

function m4(t::Triangulation; n_randpts = 100, sample_randomly = false)
	mm_dd4(t, n_randpts = n_randpts, sample_randomly = sample_randomly)
end

@testset "Markov matrix approx: version comparison" begin
	emb = invariant_gaussian_embedding(npts = 30)
	t = triang_from_embedding(Embedding(emb))

	@testset "Uniform point cloud" begin
		@time M1 = m1(t, n_randpts = 30, sample_randomly = false)
		@time M2 = m2(t, n_randpts = 30, sample_randomly = false)
		@time M3 = m2(t, n_randpts = 30, sample_randomly = false)
		@time M4 = m2(t, n_randpts = 30, sample_randomly = false)
		@time MEX = mm_p(t)

		@test all(sum(M1, 2) .≈ 1.0)
		@test all(sum(M2, 2) .≈ 1.0)
		@test all(sum(M3, 2) .≈ 1.0)
		@test all(sum(M4, 2) .≈ 1.0)
		@test all(sum(MEX, 2) .≈ 1.0)

		@show "Maximum discrepancy between exact and approx estimators:"
		@show maximum(abs(M4 .- MEX))
	end

	@testset "Random point cloud" begin
		@time M2 = mm_dd2(t, n_randpts = 10, sample_randomly = true)
		@time M4 = mm_dd4(t, n_randpts = 10, sample_randomly = true)

		@test all(sum(M2, 2) .≈ 1.0)
		@test all(sum(M4, 2) .≈ 1.0)
	end
end
