function markovmatrix(points::Array{Float64, 2},
                      image_points::Array{Float64, 2},
                      simplex_inds::Array{Int, 2})

    const n_simplices = size(simplex_inds, 1)
    const dim = size(points, 2)

    # Mismatch for the markovity of the Markov matrix is at most delta if volume_tolerance tolerance
    # is defined as follows
    const delta::Float64 = 1/10^5
    const volume_tolerance::Float64 = delta/n_simplices

    # Tolerance for similary of convex expansion coefficients of simplex vertices in simplexintersection function.
    const tolerance::Float64 = 1/10^12

    P = zeros(Float64, n_simplices, n_simplices)
    #println("There are ", n_simplices, " simplices in the triangulation.")
    for i = 1:n_simplices
        image = image_points[simplex_inds[i, :], :].'
        # Volume of the image of the simplex
        image_volume = abs(det([ones(1, dim + 1); image]))
        t1 = time_ns()
        for j = 1:n_simplices
            simplex = points[simplex_inds[j, :], :].'
            simplex_volume = abs(det([ones(1, dim + 1); simplex]))

            if (simplex_volume * image_volume > 0 && simplex_volume/image_volume > volume_tolerance)
                P[i, j] = SimplexIntersection.simplexintersection(simplex, image)/image_volume
            elseif image_volume < tolerance
                P[i, j] = 0
            end
        end
        println("Time for image simplex #", i, "\n", (time_ns() - t1)/10^9, " s")

    end

    return P
end

function markovmatrix(t::Triangulation)

    const n_simplices = size(t.simplex_inds, 1)
    const dim = size(points, 2)

    # Mismatch for the markovity of the Markov matrix is at most delta if volume_tolerance tolerance
    # is defined as follows
    const delta::Float64 = 1/10^5
    const volume_tolerance::Float64 = delta/n_simplices

    # Tolerance for similary of convex expansion coefficients of simplex vertices in simplexintersection function.
    const tolerance::Float64 = 1/10^12

    P = zeros(Float64, n_simplices, n_simplices)
    #println("There are ", n_simplices, " simplices in the triangulation.")
    for i = 1:n_simplices
        # Volume of the image of the simplex
        image_volume = t.volumes_im[i]
        for j = 1:n_simplices
            simplex_volume = t.volumes[j]
            if (simplex_volume * image_volume > 0 && simplex_volume/image_volume > volume_tolerance)
                P[i, j] = Simplices.simplexintersection(view(t.points[t.simplex_inds[j, :], :].', :, :), view(t.image_points[t.simplex_inds[i, :], :].', :, :))/image_volume
            elseif image_volume < tolerance
                P[i, j] = 0
            end
        end
        #println("Time for image simplex #", i, "\n", (time_ns() - t1)/10^9, " s")
    end

    return P
end
