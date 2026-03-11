function kinetic_monte_carlo_step!(
    szrp::SingleZRP,
    params::Parameters,
    t::Integer,
)
    for _ in 1:params.N
        i = rand(1:params.N)
        xi = szrp.particles[i]
        ni = szrp.lattice[xi]

        hop_prob = szrp.hop_rate(ni, params) * (params.dt / ni)
        if rand() < hop_prob
            hop!(szrp, params, i, xi)
        end
    end
    return nothing
end

function hop!(
    szrp::SingleZRP,
    params::Parameters,
    particle_index::Integer, 
    old_position,
)
    new_position = _new_position(old_position, szrp, params)
    szrp.particles[particle_index] = new_position
    szrp.lattice[old_position] -= 1
    szrp.lattice[new_position] += 1
    
    return nothing 
end

function _new_position(x0, szrp::SingleZRP, params::Parameters)
    x1 = x0 + rand(rng, (1, -1))
    x1 = szrp.bc(x1, params)
    return x1
end
