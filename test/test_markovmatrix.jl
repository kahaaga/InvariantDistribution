

@testset "Computing distributions" begin
    ts4::Vector{Float64} = [-1.29166, -0.0138715, -0.166332, -1.10423, -0.123669, -0.333283, -1.63198, 0.104893, -0.525271, -1.02037, -0.217897, -0.411236, -1.56985, 0.936447, 0.00540328, -0.817672, -0.66606, -0.697607, -0.848608, 0.377567, -1.11741, 1.73203, 1.36521, 0.0153917, -0.14448, -0.939471, -0.397516, -0.372424, -1.29114, 1.09698, -0.176772, 0.075443, -0.750092, -0.718689, -0.890431, -0.255283, -1.82462, 0.838895, 0.465424, -0.59368, -0.324783, -1.10423, 0.596902, -1.11463, 1.2048, 1.05063, -0.109606, -0.264874, -1.60936, -0.0519972, -0.114344, -1.08595, -0.183041, -0.424608, -1.55924, 0.751244, -0.385355, -0.199314, -1.69702, 1.39298, 1.46684, 0.00402413, -0.957294, -0.392595, -0.980378, -0.0983017, -0.65489, -1.06448, 0.360916, -1.6192, 1.50238, 1.76332, 1.29788, 0.0212719, -0.0279064, -0.396213, -1.38781, 0.0035086, -0.109723, -0.850478, -0.609453, -0.206752, -1.23392, 0.993249, -0.342366, 0.176762, -0.554082, -0.865173, -0.393957, -0.464267, -1.08186, 0.688451, -0.846321, 0.706613, -0.475809, -0.0211409, -1.70142, 1.6549, 1.78347, 0.396218, -1.10942, 1.68192, 1.43142, -0.00429295, -0.322581, -1.32556, 0.0219429, -0.0295917, -0.405318, -1.38712, 0.00763867, -0.139915, -0.948333, -0.384921, -0.386374, -1.29076, 1.1008, -0.171684, 0.0695431, -0.756431, -0.715106, -0.88513, -0.334492, -1.71579, 0.564329, -0.481763, -0.407606, -1.03805, 0.465979, -1.32223, 1.46151, 1.80423, -0.0718727, -1.85739, 1.34291, 1.99842, 0.901415, -0.51084, 0.34916, -0.688936, -0.337897, -1.23121, 0.958487, -0.403893, 0.234622, -0.620584, -0.676877, -0.742962, 0.0569321, -1.3654, 1.81456, 0.502851, -0.81959, 1.1458, 0.977908, -0.122892, -0.405326, -1.55857, 0.301308, -0.772338, -0.25191, -1.50797, 1.38587, 0.569416, -1.23556, 1.86485, 1.31952, 0.0751918, -0.774186, -0.672818, -0.980621, -0.165931, -1.87667, 0.260215, -0.831583, -0.213941, -1.70374, 1.35865, 1.42889, 0.0350163, -0.938797, -0.402145, -1.17729, -0.00212732, -0.18619, -1.11735, -0.092832, -0.303358, -1.56958, -0.028755, -0.11704, -0.996922, -0.315869, -0.452005, -1.30119, 1.12457, -0.135106, 0.0117999, -0.924947, -0.447779, -0.992671]
    SS = markovmatrix(ts4, 3, 2)

    @testset "Markov matrix" begin
        @test all(sum(SS.P, 2) .≈ 1.0)
    end

    invmeasure, inds_nonzero_simplices = invariantdist(SS.P)

    @testset "Invariant measure" begin
        # There must be at least one simplex with nonzero measure; otherwise, the invariant
        # distribution cannot sum to 1.
        @test size(inds_nonzero_simplices, 1) >= 1

        # The invariant distribution must sum to 1 within reasonable tolerance.
        @test sum(invmeasure) ≈ 1.0 atol = 1/10^6
    end


    @testset "Find nonempty bins" begin
        nonzero_simplices = SS.triangulation_indices[inds_nonzero_simplices, :]
        vertices = SS.triangulation_vertices
        centroids, radii = SimplexSplitting.centroids_radii2(vertices, nonzero_simplices)
        #@show centroids
        @testset "n_bins = $n" for n in 15:1:15
            nonempty_bins, measure = get_nonempty_bins(centroids, invmeasure, [n, n, n])
            @test sum(measure) ≈ 1.0 atol = 1/10^6
            @test !any(nonempty_bins .== 0)

            @testset "joint" begin
                joint = jointdist(nonempty_bins, measure)
                #@show joint
                @show size(joint)
            end

            @testset "XY" begin
                #Py, Pxy, Pxz, JY, JXY, JXZ = marginaldists(nonempty_bins, measure)
                Py, Pxy, Pxz = marginaldists(nonempty_bins, measure)
                @show size(Py)
                @show size(Pxy)
                @show size(Pxz)
            end
        end
    end



end
