using InvariantDistribution, SimplexSplitting, StaticArrays, Distributions,
	BenchmarkTools

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

@testset "mm_discrete" begin
	t = triang_from_embedding(Embedding(invariant_gaussian_embedding(npts = 8, covariance = 0.5, tau = 1)))
	println("\nTriangulation has $(size(t.simplex_inds, 1)) simplices.")
	npts = 100
	println("\tApproximating M by uniform sampling of at least $(npts) points.")
	print("\t"); @time M1 = m1(t, n_pts = npts, sample_randomly = false)
	print("\t"); @time M2 = m2(t, n_pts = npts, sample_randomly = false)
	print("\t"); @time M3 = m3(t, n_pts = npts, sample_randomly = false)

	t = triang_from_embedding(Embedding(invariant_gaussian_embedding(npts = 30, covariance = 0.5, tau = 1)))
	println("\nTriangulation has $(size(t.simplex_inds, 1)) simplices.")
	println("\tApproximating M by uniform sampling of at least $(npts) points.")
	print("\t"); @time M1 = m1(t, n_pts = npts, sample_randomly = false)
	print("\t"); @time M2 = m2(t, n_pts = npts, sample_randomly = false)
	print("\t"); @time M3 = m3(t, n_pts = npts, sample_randomly = false)

	@test all(sum(M1, 2) .≈ 1)
	@test all(sum(M2, 2) .≈ 1)
	@test all(sum(M3, 2) .≈ 1)
	@test M1 == M2 == M3

	println("\n\tApproximating M by random sampling of exactly $(npts) points.")
	print("\t"); @time M1 = m1(t, n_pts = npts, sample_randomly = true)
	print("\t"); @time M2 = m2(t, n_pts = npts, sample_randomly = true)
	print("\t"); @time M3 = m3(t, n_pts = npts, sample_randomly = true)
	@test all(sum(M1, 2) .≈ 1)
	@test all(sum(M2, 2) .≈ 1)
	@test all(sum(M3, 2) .≈ 1)
end
