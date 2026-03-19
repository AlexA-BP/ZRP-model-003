function naive_monte_carlo_step!(
    state::NSpeciesState,
    model::ModelParameters,
    modelfunc::ModelFunctions,
    t::Integer,
    rng::AbstractRNG,
)
    for _ in 1:model.N
        i = rand(rng, 1:model.N)
        xi = state.particles[i]
        si = state.species[i]

        hop_prob = get_hop_prob(xi, si, state, model, modelfunc) 
        if rand(rng) < hop_prob
            hop!(state, model, modelfunc, i, si, xi, rng)
        end
    end
    return nothing
end

function get_hop_prob(
    xi::Integer, 
    si::Integer,
    state::NSpeciesState,
    model::ModelParameters, 
    modelfunc::ModelFunctions,
)
    n_A_i = state.lattice[xi, 1]
    n_B_i = state.lattice[xi, 2]
    
    single_hop_rate = modelfunc.hop_rate(si, n_A_i, n_B_i, model.alpha, model.chi)

    return state.lattice[xi, si] * single_hop_rate * model.dt
end

function hop!(
    state::NSpeciesState,
    model::ModelParameters,
    modelfunc::ModelFunctions,
    particle_index::Integer, 
    particle_species::Integer,
    old_position::Real,
    rng::AbstractRNG,
)
    new_position = _new_position(old_position, modelfunc.bc, model, rng)
    state.particles[particle_index] = new_position
    state.lattice[old_position, particle_species] -= 1
    state.lattice[new_position, particle_species] += 1
   
    return nothing 
end

function _new_position(x0, bc::Function, model::ModelParameters, rng::AbstractRNG)
    x1 = x0 + rand(rng, model.hop_directions)
    x1 = bc(x0, model)
    return x1
end
