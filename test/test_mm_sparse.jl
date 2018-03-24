
@testset "Invariant distribution tests" begin
   @testset "Embedding dim=$E" for E in 3:3
      @testset "Time series length: $ts_length" for ts_length in [10]
         @testset "Rep #$rep" for rep in 1:1
            println("E=",E, "\tlength(ts)=", ts_length, "\n")
            t_start = time_ns()

            # Create an invariant embedding
            embedding = gaussian_embedding(ts_length + E + 10)

            t = triang_from_embedding(embedding)

            # Split the triangulation such that all simplices have radii below some fraction
            # of the largest simplex radius in the triangulation.
            #print("Splitting triangulation...\n")
            #print("Initial number of simplices = ", size(t.simplex_inds, 1), "\n")


            #SimplexSplitting.refine_variable_k!(t, maximum(t.radii) - (maximum(t.radii) - mean(t.radii))/2)
            #print("Final number of simplices = ", size(t.simplex_inds, 1), "\n")

            #print("The number of simplices increased from...\n")

            # Create Markov matrix.
            #print("Testing sparse markov matrix computation\n")
            #t_sparse = time_ns()
            #M = mm_sparse(t)
            #print("mm_sparse took ", (time_ns() - t_sparse)/10^9, " seconds\n\n")
            print("Testing parallel markov matrix computation\n")
            t_parallel = time_ns()
            Mp = mm_p(t)
            print("mm_p took ", (time_ns() - t_parallel)/10^9, " seconds\n\n")

            t_markov = time_ns()

            # The column sums of the Markov matrix all must be 1. Otherwise, we haven't
            # accounted for the entire volume of the triangulation.
            @testset "Markov matrix" begin
               println(sum(Mp, 2))
               @test all(sum(Mp, 2) .≈ 1.0)
            end


            invmeasure, inds_nonzero_simplices = invariantdist(Array(Mp))
            t_measure = time_ns()
            #println("E=",E, "\tlength(ts)=", ts_length, "\trep=", rep, " | Markov: ", (t_markov - t_start)/10^9, " seconds for ", size(t.simplex_inds, 1), " simplices with ", nnz(M), " nonzero intersections, invmeasure: ", (t_measure - t_markov)/10^9, " seconds")

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
