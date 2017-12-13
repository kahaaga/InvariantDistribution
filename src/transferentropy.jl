"""
Compute transfer entropy from joint and marginal distributions
"""
function te(nonempty_bins::Array{Int, 2}, joint::Vector{Float64},
    marginals::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}})

    TE = 0

    JY = indexin()

    for i = 1:size(nonempty_bins, 1)
        Pxyz = joint[i]
        Py =
    end

    Pxyz=PJoint(i);
           Py=PY(JY(i));
           Pxy=PXY(JXY(i));
           Pyz=PYZ(JYZ(i));
           te=te+Pxyz*log(Pxyz*Py/(Pyz*Pxy));
end
