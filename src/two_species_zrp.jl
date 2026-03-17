#=
Czech-list
    [ ] basic structure
    [ ] file handling
        [ ]
        [ ]
    [ ] simulation  
=#
using Random

include("./sim/TwoSpeciesZRP/Structs.jl")
include("./sim/HoppingRates.jl")
include("./sim/BoundaryConditions.jl")
include("./sim/filehandling.jl")
include("./sim/kinetic_monte_carlo.jl")

import .HoppingRates: hop_rates
using .BoundaryConditions


function main_dev(;
    # function paramters, for now hardcoded
    num_particles_A=5,
    num_particles_B=6,
    system_size=10,
    bc="p",
    alpha=1., 
    chi=3., 
    tot_timesteps=100,
    chunk_size=10,
    time_steps_between_snapshots=5,
    fname="./data/two_species.h5"
)


    # const rng = Xoshiro(2) 
    rng = Random.default_rng()
    this_hop_rate, get_dt = hop_rates["simple_weak_strong_nr_hop_rate"]


    # handle function input and create necessary parameters for everything
    model = ModelParameters(
        [num_particles_A, num_particles_B],
        2,
        num_particles_A + num_particles_B,
        system_size,
        tot_timesteps,
        get_dt(num_particles_A, num_particles_B, alpha, chi),
        bc,
        alpha,
        chi,
        (1, -1),
    )
    state = NSpeciesState(model, rng)
    modelfunc = ModelFunctions(this_hop_rate, bcs[bc])
    chunk = Chunks(chunk_size, time_steps_between_snapshots, tot_timesteps)    

    # setup hdf5
    setup_hdf5(fname, model, state, chunk)
    
    run_and_write_chunked_simulation(
        fname,
        model,
        state,
        modelfunc,
        chunk,
        rng,
    )
    return (model, state)

end
