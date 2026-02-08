module BoundaryConditions

export bcs

bcs = Dict()

function bc_periodic(x, params)
    if x < 1
        return params.L
    elseif x > params.L
        return one(x)
    end
    return x
end
bcs["p"] = bc_periodic
bcs["periodic"] = bc_periodic


function bc_reflective(x, params)
    if x < 1
        return 2*one(x)
    elseif x > params.L
        return params.L - one(x)
    end
    return x
end
bcs["r"] = bc_reflective
bcs["reflective"] = bc_reflective

end