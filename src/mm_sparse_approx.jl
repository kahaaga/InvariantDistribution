function mm_sparse(t::SimplexSplitting.Triangulation)

    const n_simplices = size(t.simplex_inds, 1)
    const dim = size(t.points, 2)

    # Mismatch for the markovity of the Markov matrix is at most delta if volume_tolerance
    # tolerance as below
    const delta::Float64 = 1/10^5
    const voltol::Float64 = delta/n_simplices

    # Tolerance for similary of convex expansion coefficients of simplex vertices in simplexintersection function.
    const convex_params_tol::Float64 = 1/10^12

    Is = Vector{Int}(0)
    Js = Vector{Int}(0)
    intvols = Vector{Float64}(0) # intersecting volumes

    kdtree = KDTree(t.centroids)

    for i in 1:n_simplices
        imvol = t.volumes_im[i]
        t1 = time_ns()
        # Consider only simplices in the original space that are close to any of the
        # vertices or the centroid of this image simplex
        inds_to_check =

        @elapsed for j in 1:n_simplices
            vol = t.volumes[j]

            if vol * imvol > 0 && (vol/imvol) > voltol
                intvol = Simplices.simplexintersection(t.points[t.simplex_inds[j, :], :].', t.impoints[t.simplex_inds[i, :], :].')/imvol

                if intvol > 0
                    push!(Is, i)
                    push!(Js, j)
                    push!(intvols, intvol)
                # else
                #     if i == n_simplices && j == n_simplices
                #         push!(Is, i)
                #         push!(Js, j)
                #         push!(intvols, 0)
                #     end
                end
            end
        end
        println("Image simplex #",i, ":", (time_ns() - t1)/10^9, " s elapsed")

    end
    return sparse(Is, Js, intvols, n_simplices, n_simplices)
end
