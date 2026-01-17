#=
    A bunch of small utilities.
=#



function assign_hyperslab!(A, B, inds)
    selectdim(A, ndims(A), inds) .= B
    return nothing
end

function assign_hyperslab!(A::Matrix, v::Vector, inds)
    A[:, inds] = v
    return nothing
end