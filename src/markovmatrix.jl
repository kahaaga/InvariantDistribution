function markovmatrix(points, image_points, simplex_inds)

    const n_simplices = size(simplex_inds, 1)
    const dim = size(points, 2)

    # Mismatch for the markovity of the Markov matrix is at most delta if volume_tolerance tolerance
    # is defined as follows
    const delta::Float64 = 1/10^5
    const volume_tolerance::Float64 = delta/n_simplices

    # Tolerance for similary of convex expansion coefficients of simplex vertices in simplexintersection function.
    const tolerance::Float64 = 1/10^12

    P = zeros(n_simplices, n_simplices)

    for i = 1:n_simplices
        image = image_points[simplex_inds[i, :], :].'

        # Volume of the image of the simplex
        image_volume = abs(det([ones(1, dim + 1); image]))

        # test if
        if image_volume < tolerance
            # Loop over image vertices
            for col in 1:size(image, 2)

              # Vertices in triangulation
                for vertex in unique_vertex_inds
                    v = triangulation[1][vertex, :]
                    if isapprox(image[:, col], v)
                        println("image simplex #", i, " = ", image, " (image_vol ",
                              image_volume, ")\n shares vertex with point #", vertex,
                              " in triangulation (coordinates = ", triangulation[1][vertex, :], ")\n")
                    end
                end
            end
        end

        for j = 1:n_simplices
            simplex = points[simplex_inds[j, :], :].'
            simplex_volume = abs(det([ones(1, dim + 1); simplex]))

            if (simplex_volume * image_volume > 0 && simplex_volume/image_volume > volume_tolerance)
                P[i,j] = SimplexIntersection.simplexintersection(simplex, image)/image_volume
            end
        end
    end

    return P
end
function markovmatrixp(points, image_points, simplex_inds)

    const n_simplices = size(simplex_inds, 1)
    const dim = size(points, 2)

    # Mismatch for the markovity of the Markov matrix is at most delta if volume_tolerance tolerance
    # is defined as follows
    const delta::Float64 = 1/10^5
    const volume_tolerance::Float64 = delta/n_simplices

    # Tolerance for similary of convex expansion coefficients of simplex vertices in simplexintersection function.
    const tolerance::Float64 = 1/10^12

    P = SharedArray{Float64}(n_simplices, n_simplices)

    @parallel for i = 1:n_simplices
        image = image_points[simplex_inds[i, :], :].'

        # Volume of the image of the simplex
        image_volume = abs(det([ones(1, dim + 1); image]))

        # test if
        if image_volume < tolerance
            # Loop over image vertices
            for col in 1:size(image, 2)

              # Vertices in triangulation
                for vertex in unique_vertex_inds
                    v = triangulation[1][vertex, :]
                    if isapprox(image[:, col], v)
                        println("image simplex #", i, " = ", image, " (image_vol ",
                              image_volume, ")\n shares vertex with point #", vertex,
                              " in triangulation (coordinates = ", triangulation[1][vertex, :], ")\n")
                    end
                end
            end
        end

        for j = 1:n_simplices
            simplex = points[simplex_inds[j, :], :].'
            simplex_volume = abs(det([ones(1, dim + 1); simplex]))

            if (simplex_volume * image_volume > 0 && simplex_volume/image_volume > volume_tolerance)
                P[i,j] = SimplexIntersection.simplexintersection(simplex, image)/image_volume
            end
        end
    end

    return fetch(P)
end
