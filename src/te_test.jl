
"""
Wrapper allowing the computation of transfer entropy directly from a time series, given an
embedding dimension and an embedding lag. Only for testing. For real use cases, use the version
of the function taking an embedding.
"""
function te_test(ts::Array{Float64}, embeddim, embedlag, binsizes)
    embedding = embed(ts, embeddim, embedlag)

    if invariantset(embedding)
      println("The embedding forms an invariant set. Continuing.")
    else
      println("The embedding does not form an invariant set. Quitting.")
      return
    end

    # Embed using all points except the last (to allow projection of all vertices
    # in the triangulation).
    points, simplex_inds = triangulate(embedding[1:end-1, :])
    image_points = embedding[2:end, :]

    P = markovmatrix(points, image_points, simplex_inds)

    invmeasure, inds_nonzero_simplices = invariantdist(P)


    nonzero_simplices = simplex_inds[inds_nonzero_simplices, :]
    vertices = points
    centroids, radii = SimplexSplitting.centroids_radii2(vertices, nonzero_simplices)

    TE = zeros(Float64, length(binsizes))

    for i = 1:length(binsizes)
        n = binsizes[i]
        # Find bins with nonzero measure.
        nonempty_bins, measure = get_nonempty_bins(centroids, invmeasure, [n, n, n])
        unique_nonempty_bins = unique(nonempty_bins, 1)

        # Compute joint and marginal distributions.
        joint = jointdist(nonempty_bins, measure)
        Py, Pxy, Pxz = marginaldists(unique_nonempty_bins, measure)
        marginals = Py, Pxy, Pxz

        # Compute transferentropy
        TE[i] = transferentropy(unique_nonempty_bins, joint, marginals)
    end

    return TE
end
