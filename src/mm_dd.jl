using Simplices, SimplexSplitting, StaticArrays

function potentially_intersecting_simplices(t::Triangulation, image_i::Int)
    inds_potential_simplices = Int[]

    n_simplices = length(t.radii)

    @inbounds for i = 1:n_simplices
        dist_difference = ((t.centroids_im[image_i] - t.centroids[i]).' *
                            (t.centroids_im[image_i] - t.centroids[i]) - (t.radii_im[image_i] + t.radii[i])^2)[1]
        if dist_difference < 0
            push!(inds_potential_simplices, i)
        end
    end
    return inds_potential_simplices
end


function contains_point!(signs, s_arr, s, point)

    # Keep track of signs of determinants. If all have the same sign, the point
    # lies inside the simplex. If not, it lies outside.
    n_pts = size(s, 1)
    # Compute the first sign
    @views s_arr[1, 1:(end-1)] = point
    @views s_arr[2:end, 1:(end-1)] = s[2:end, :]

    signs[1] = sign(det(s_arr))

    for i = 2:(n_pts - 1) # Check remaining signs and stop if sign changes
        @views s_arr[1:(i - 1), 1:(end-1)] = s[1:(i - 1), :]
        @views s_arr[i, 1:(end-1)] = point
        @views s_arr[(i + 1):end, 1:(end-1)] = s[(i + 1):end, :]

        signs[i] = sign(det(s_arr))
        if !(signs[i] == signs[i - 1])
           return false
        end
    end

    @views s_arr[1:(end - 1), 1:(end-1)] = s[1:(end - 1), :]
    @views s_arr[end, 1:(end-1)] = point

    signs[end] = sign(det(s_arr))
    if !(signs[end] == signs[end - 1])
       return false
    end
    return true
end


function mm_dd(t::Triangulation;
    n_randpts::Int = 100,
    dist::Distributions.Distribution = Distributions.Uniform(0, 1),
    sample_randomly::Bool = false)

    n_simplices = size(t.simplex_inds, 1)
    n_vertices = size(t.simplex_inds, 2)
    dim = n_vertices - 1
    simplices = Size(n_simplices, dim + 1, dim)(zeros(n_simplices, dim + 1, dim))
    imsimplices = Size(n_simplices, dim + 1, dim)(zeros(n_simplices, dim + 1, dim))

    for i in 1:n_simplices
        simplices[i, :, :] = t.points[t.simplex_inds[i, :], :]
        imsimplices[i, :, :] = t.impoints[t.simplex_inds[i, :], :]
    end

    convex_coeffs = get_coeffs(dim, n_randpts, sample_randomly).'
    n_coeffs = maximum(size(convex_coeffs))

    pt = Size(dim)(zeros(Float64, dim))
    cache_arr = Size(dim + 1, dim + 1)(zeros(Float64, 4, 4))
    cache_arr[:, dim+1] = ones(dim+1)

    s_arr = Size(dim, dim + 1)(zeros(Float64, dim, dim + 1))
    signs = Size(dim + 1)(zeros(Float64, dim + 1))

    M = zeros(Float64, n_simplices, n_simplices)

    for i in 1:n_simplices
        inds = potentially_intersecting_simplices(t, i)
        image = imsimplices[i, :, :]

        for c in 1:n_coeffs
            @views pt = convex_coeffs[c, :].' * image
            for j in inds
                simplex = @views simplices[j, :, :]
                if contains_point!(signs, cache_arr, simplex, pt)
                    M[j, i] += 1.0
                end
            end
        end
    end
    return M.' ./ n_coeffs
end
