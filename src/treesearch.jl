using NearestNeighbors

dim = 3
data = rand(3, 10^4)
k = dim * dim
points = rand(3, 4)

kdtree = KDTree(data)
idxs, dists = knn(kdtree, points, k, true)


function tree_example(;n_simplices::Int = 1000, dim::Int = 3, k::Int = dim + 1)
    # Fake a triangulation by simulation simplex centroids
    centroids = rand(dim, n_simplices)
    points = rand(dim, (dim + 1) + 10)

    kdtree = KDTree(centroids)
    idxs, dists = knn(kdtree, points, k, true)
    return idxs, dists
end



@time tree_example(3000, 3, 4)
