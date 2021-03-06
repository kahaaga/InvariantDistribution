
function mm_dd3(t::Triangulation, simplices, convex_coeffs)

    # Some constants used throughout the funciton
    n_simplices = size(t.simplex_inds, 1)
    dim = size(t.points, 2)

    #=
    # update number of points in case a regular grid of points was employed
    # (if so, because the number of subsimplices generated by the
    # shape-preserving splitting depends only on the dimension of the space,
    # there will be more points than we asked for).
    =#
    n_coeffs = maximum(size(convex_coeffs))


    # Pre-allocated arrays (SizedArrays, for efficiency)
    pt          = Size(dim)(zeros(Float64, dim))
    s_arr       = Size((dim+1)^2)(zeros(Float64, (dim+1)^2))
    signs       = Size(dim + 1)(zeros(Float64, dim + 1))

    # Re-arrange simplices so that look-up is a bit more efficient
    simplex_arrs = Vector{Array{Float64, 2}}(n_simplices)
    imsimplex_arrs = Vector{Array{Float64, 2}}(n_simplices)
    for i in 1:n_simplices
        simplex_arrs[i] = t.points[t.simplex_inds[i, :], :]
        imsimplex_arrs[i] = t.impoints[t.simplex_inds[i, :], :]
    end

    # The Markov matrix
    M = zeros(Float64, n_simplices, n_simplices)

    for i in 1:n_simplices
        inds = potentially_intersecting_simplices2(t, i)
        is = imsimplex_arrs[i]

        @views  for k in 1:n_coeffs
            pt = convex_coeffs[:, k].' * is

            #@show i, k, pt
            for j in inds
                sx = simplices[:, j]
                if contains_point_lessalloc!(signs, s_arr, sx, pt, dim)
                    M[j, i] += 1.0
                end
            end

        end
    end
    return M.' ./ n_coeffs
end
