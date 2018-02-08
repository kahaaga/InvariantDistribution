
@testset "Invariant distribution tests" begin
   @testset "Embedding dim=$E" for E in 3:5
      @testset "Time series length: $ts_length" for ts_length in 50:50:50
         @testset "Rep #$rep" for rep in 1:1
            t_start = time_ns()

            dist = Normal()
            # Create an invariant embedding
            embedding = invariant_embedding(; dist = dist, npts = ts_length + E - 1)

            t = triang_from_embedding(embedding)

            # Split the triangulation such that all simplices have radii below the mean
            # radius of the simplices.
            SimplexSplitting.refine_variable_k!(t, mean(t.radii))

            # Create Markov matrix.
            M = mm_sparse(t)
            t_markov = time_ns()

            # The column sums of the Markov matrix all must be 1. Otherwise, we haven't
            # accounted for the entire volume of the triangulation.
            @testset "Markov matrix" begin
               @test all(sum(M, 2) .≈ 1.0)
            end


            invmeasure, inds_nonzero_simplices = invariantdist(M)
            t_measure = time_ns()
            println("E=",E, "\tlength(ts)=", ts_length, "\trep=", rep, " | Markov: ", (t_markov - t_start)/10^9, " seconds for ", size(t.simplex_inds, 1), " simplices with ", nnz(M), " nonzero intersections, invmeasure: ", (t_measure - t_markov)/10^9, " seconds")

            # There must be at least one simplex with nonzero measure; otherwise, the
            # invariant distribution cannot sum to 1. Also, the invariant distribution must
            # sum to 1 within reasonable tolerance.
            @testset "Invariant measure" begin
               @test size(inds_nonzero_simplices, 1) >= 1
               @test sum(invmeasure) ≈ 1.0 atol = 1/10^6
            end
         end
      end
   end
end
