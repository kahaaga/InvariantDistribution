
"""
Return a triangulation of the canonical simplex. Used for all the invariant
maps in this file.
"""
function canonical_simplex_triangulation(; dim::Int = 3, split_factor::Int = 3)
    # Define vertices of canonical simplex
    canonical_simplex_vertices = zeros(dim + 1, dim)
    canonical_simplex_vertices[2:(dim+1), :] = eye(dim)
    simplex_indices = zeros(Int, 1, dim + 1)
    simplex_indices[1, :] = collect(1:dim+1)

    refined = refine_triangulation(
        canonical_simplex_vertices,
        simplex_indices, [1],
        split_factor)

    points, simplex_inds = refined[1], refined[2]
    centroids, radii = SimplexSplitting.centroids_radii2(points, simplex_inds)
    orientations = SimplexSplitting.orientations(points, simplex_inds)
    volumes = abs.(orientations)

    Triangulation(
        points = points,
        simplex_inds = simplex_inds,
        centroids = centroids,
        radii = radii,
        orientations = orientations,
        volumes = volumes
    )
end


"""
    invariant_map_bull(t::Triangulation)

Invariant map given by

Ψ(x, y, z) = (x/2 + z/2 * (1-y), (y+x)/2, 1-z),

where x, y and z are columns of the t.points array. These points are guaranteed to lie in the unit cube, because each component <= 1.
"""
function invariant_map_bull(t::Triangulation)
    hcat(
        (t.points[:, 1]/2 + t.points[:, 3]/3) .* (1 - t.points[:, 2]),
        (t.points[:, 1] + t.points[:, 2])/2,
        1 - t.points[:, 3]
    )
end



function invariant_bullmap_on_cube(;dim::Int = 3, split_factor::Int = 4)

    if dim == 3
        t = canonical_simplex_triangulation(dim = dim, split_factor = split_factor)
        t.impoints = invariant_map_bull(t)
        t.centroids_im, t.radii_im = SimplexSplitting.centroids_radii2(t.impoints, t.simplex_inds)
        t.orientations_im = SimplexSplitting.orientations(t.impoints, t.simplex_inds)
    else
        error("map not defined for dim = ", dim)
    end

    return t
end
