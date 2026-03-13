function kinetic_monte_carlo_step!(
    state::TwoSpeciesState,
    prm::Parameters,
    modelfunc::ModelFunctions,
    t::Integer,
)
    for _ in 1:prm.N
        i = rand(1:prm.N)
        if i <= prm.N_A
            species = "A"
            xi = state.particles_A[i]
            n_S = state.lattice_A[xi]
        elseif i > prm.N_A
            species = "B"
            xi = state.particles_B[i-prm.N_A]
            n_S = state.lattice_B[xi]
        end
        n_A = state.lattice_A[xi]
        n_B = state.lattice_B[xi]

        hop_prob = n_S*modelfunc.hop_rate(species, n_A, n_B, prm.alpha, prm.chi) * prm.dt
        if rand() < hop_prob
            hop!(state, prm, modelfunc, i, xi)
        end
    end
    return nothing
end

function hop!(
    state::TwoSpeciesState,
    prm::Parameters,
    modelfunc::ModelFunctions,
    particle_index::Integer, 
    old_position::Real,
)
    new_position = _new_position(old_position, modelfunc.bc, prm)
    if particle_index <= prm.N_A
        state.particles_A[particle_index] = new_position
        state.lattice_A[old_position] -= 1
        state.lattice_A[new_position] += 1
    elseif particle_index > prm.N_B
        state.particles_B[particle_index-prm.N_A] = new_position
        state.lattice_B[old_position] -= 1
        state.lattice_B[new_position] += 1
    end
    
    return nothing 
end

function _new_position(x0, bc::Function, prm::Parameters)
    x1 = x0 + rand(rng, (1, -1))
    x1 = bc(x1, prm)
    return x1
end
