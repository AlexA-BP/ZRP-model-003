#=
    A bunch of small utilities.
=#



function assign_hyperslab!(A, B, inds)
    selectdim(A, ndims(A), inds) .= B
    return nothing
end