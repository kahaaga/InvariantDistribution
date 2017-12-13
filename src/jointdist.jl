"""
Compute the joint distribution of variables from a triangulation of an embedding
of the variables X and Y {(x(t), x(t-tau), y(t))}.

`centroids::Array{Float, 2}` are centroids of the simplices forming the triangulation; this
is an with size n_simplices x dim.
`invariantdist::Array{Float64, 1}` is the invariant distribution on the triangulation
`δ::Vector{Float64}` are the number of bins along each dimension. The algorithm
    divides each axis into equidistant bins from min(centroids[:, i]) to max(centroids[:, i])
"""
function get_nonempty_bins(centroids::Array{Float64, 2},
                    invariantdist::Array{Float64, 1},
                    δ::Vector{Int})

    n_simplices = size(centroids, 1)
    dim = size(centroids, 2)

    # The number of grids specified must match the number of dimensions.
    @assert dim == length(δ)

    # Initialise matrix holding the non-zero entries of the joint distribution.
    # Each row contains the indices
    nonempty_bins = zeros(Int, n_simplices, dim)

    # Initialise column vector holding the invariant densities of the boxes corresponding
    # to rows in 'joint'.
    densities = zeros(Float64, n_simplices)

    # For each simplex, find the index of the bin its centroid lies in for each of the
    # dimensions. For this, we need to know where to locate the bins.
    startvals = [minimum(centroids[:, i]) for i in 1:dim]
    endvals = [maximum(centroids[:, i]) for i in 1:dim]

    ranges = [(abs(startvals[i]) + abs(endvals[i])) for i in 1:dim]

    for i = 1:n_simplices
        densities[i] = invariantdist[i]
        for j = 1:dim
            stepsize = ranges[j] / δ[j]
            pos_along_range = abs(centroids[i, j] - startvals[j])
            nonempty_bins[i, j] = floor(Int, pos_along_range / stepsize) + 1
        end
    end

    return nonempty_bins, densities
end


function jointdist(nonempty_bins::Array{Int, 2}, densities::Vector{Float64})

    unique_bins = unique(nonempty_bins, 1)
    Pjoint = zeros(Float64, size(unique_bins))

    for i = 1:size(unique_bins, 1)
        inds = find(all(nonempty_bins .== unique_bins[i, :].', 2))
        Pjoint[i] = sum(densities[inds])
    end

    return Pjoint

end


function marginaldists(nonempty_bins::Array{Int, 2}, densities::Vector{Float64})
    dim = size(nonempty_bins, 2)


    # PY

    # For each unique X2, sum the measure of all bins which has that X2 coordinate
    X2s = nonempty_bins[:, 2]
    X2s = reshape(X2s, length(X2s), 1)
    X1X2s = nonempty_bins[:, 1:2]
    X2X3s = nonempty_bins[:, 2:3]
    unique_X2s = unique(X2s, 1)
    unique_X1X2s = unique(X1X2s, 1)
    unique_X2X3s = unique(X2X3s, 1)
    JX2 = indexin(X2s, unique_X2s)
    JX1X2 = indexin_rows(X1X2s, unique_X1X2s)
    JX2X3 = indexin_rows(X2X3s, unique_X2X3s)


    PX2 = zeros(Float64, size(unique_X2s, 1))
    for i = 1:size(unique_X2s, 1)
        inds = find(JX2 .== i)
        PX2[i] = sum(densities[inds])
    end
    Px2 = zeros(Float64, size(X2s, 1))

    for i = 1:size(X2s, 1)
        Px2[i] = PX2[JX2[i]]
    end

    PX1X2 = zeros(Float64, size(unique_X1X2s, 1))
    for i = 1:size(unique_X1X2s, 1)
        inds = find(JX1X2 .== i)
        PX1X2[i] = sum(densities[inds])
    end
    Px1x2 = zeros(Float64, size(X1X2s, 1))

    for i = 1:size(X1X2s, 1)
        Px1x2[i] = PX1X2[JX1X2[i]]
    end


    PX2X3 = zeros(Float64, size(unique_X2X3s, 1))
    for i = 1:size(unique_X2X3s, 1)
        inds = find(JX2X3 .== i)
        @show inds
        PX2X3[i] = sum(densities[inds])
    end
    Px2x3 = zeros(Float64, size(X2X3s, 1))

    for i = 1:size(X2X3s, 1)
        Px2x3[i] = PX2X3[JX2X3[i]]
    end

    return PX2, PX1X2, PX2X3

end
