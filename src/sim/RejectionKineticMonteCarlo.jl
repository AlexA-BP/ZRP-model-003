function rejection_kinetic_monte_carlo_step!(
    state::NSpeciesState,
    model::ModelParameters,
    modelfunc::ModelFunctions,
    rng::AbstractRNG,
)
    i = rand(rng, 1:model.N)
    # check if move gets accepted
    if rand(rng) < state.trans_rates[i]/state.r0


        old_x = state.particles[i]
        new_x = _new_position(old_x, modelfunc.bc, model, rng)
        spec = state.species[i]

        # update particle position
        hop!(state, i, old_x, new_x)

        # update the transition rates of the two sites, which changed
        state.trans_rates[old_x, spec] = modelfunc.hop_rate(
            state, model, old_x, spec
        )
        state.trans_rates[new_x, spec] = modelfunc.hop_rate(
            state, model, new_x, spec
        )

        # since r0 might've changed, potentially update
        state.r0 = maximum(state.trans_rates)

        println(state.r0)
        # determine time step and update current time
        state.dt = 1/(2*model.N * state.r0) * log(1/rand(rng))
        state.t += state.dt

    end

    return nothing
end

function hop!(
    state::NSpeciesState,
    index,
    old_x,
    new_x,
)
    spec = state.species[index]
    state.particles[index] = new_x
    state.lattice[old_x, spec] -= oneunit(old_x)
    state.lattice[new_x, spec] += oneunit(old_x)

    return nothing
end

function _new_position(
    x0::Integer, bc::Function, model::ModelParameters, rng::AbstractRNG,
)
    x1 = x0 + rand(rng, model.hop_directions)
    x1 = bc(x0, model)
    return x1
end