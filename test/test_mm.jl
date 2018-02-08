
tic()
@testset "Invariant distribution tests" begin
   dist = Normal()
   @testset "Time series length: $ts_length" for ts_length in 50:50:300
      @testset "Rep #$rep" for rep in 1:3

         is_invariant = false
         TS = Float64[]
         embedding = Array{Float64}(0, 0)
         while (!is_invariant)
            ts = rand(dist, ts_length)
            e = embed(ts, 3, 2)
            if invariantset(e)
               is_invariant = true
               TS = ts
               embedding = e
            end
         end

         println("length(ts) = ", ts_length, "\trep=", rep, "\tThe embedding forms an invariant set")
         isinvariant = true

         # Embed using all points except the last (to allow projection of all vertices
         # in the triangulation).
         points, simplex_inds = triangulate(embedding[1:end-1, :])
         image_points = embedding[2:end, :]

         @testset "Simplex vertices matches embedding[1:end-1, :]" begin
            @test points == embedding[1:end-1, :]
         end

         P = markovmatrix(points, image_points, simplex_inds)

         @testset "Markov matrix" begin
            @test all(sum(P, 2) .≈ 1.0)
         end

         invmeasure, inds_nonzero_simplices = invariantdist(P)


         @testset "Invariant measure" begin
            # There must be at least one simplex with nonzero measure; otherwise, the invariant
            # distribution cannot sum to 1.
            @test size(inds_nonzero_simplices, 1) >= 1

            # The invariant distribution must sum to 1 within reasonable tolerance.
            @test sum(invmeasure) ≈ 1.0 atol = 1/10^6
         end
            #
            # nonzero_simplices = simplex_inds[inds_nonzero_simplices, :]
            # vertices = points
            # centroids, radii = SimplexSplitting.centroids_radii2(vertices, nonzero_simplices)
            # #@show size(centroids), size(invmeasure)
            # n = 3
            # nonempty_bins, measure = get_nonempty_bins(centroids, invmeasure, [n, n, n])
            #
            # unique_nonempty_bins = unique(nonempty_bins, 1)
            #
            # @testset "Bin measure is 1" begin
            #    @test sum(measure) ≈ 1.0 atol = 1/10^6
            # end
            #
            # @testset "All centroids are assigned to bins" begin
            #    @test !any(nonempty_bins .== 0)
            # end
            #
            #
            # @testset "Joint distribution" begin
            #    joint = jointdist(nonempty_bins, measure)
            #    @show joint
            # end
            #
            # @testset "Marginal distribution" begin
            #    Py, Pxy, Pxz = marginaldists(unique_nonempty_bins, measure)
            #    @show Py
            #    @show Pxy
            #    @show Pxz
            # end
            #
            # joint = jointdist(nonempty_bins, measure)
            # Py, Pxy, Pxz = marginaldists(unique_nonempty_bins, measure)
            #
            # #@show size(joint), size(Py), size(Pxy), size(Pxz)
            #
            # @testset "Joint and marginals are the same size" begin
            #    sizes = [size(distribution, 1) for distribution in (joint, Py, Pxy, Pxz)]
            #    @test length(unique(sizes)) == 1
            # end
      end
   end
end
toc()
