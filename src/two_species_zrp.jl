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
    num_particles_A=5000,
    num_particles_B=5000,
    system_size=10,
    bc="p",
    alpha=0., 
    chi=3., 
    phys_time=10.,
    chunk_size=10,
    time_steps_between_snapshots=1,
    fname="./data/two_species.h5"
)


    # const rng = Xoshiro(2) 
    rng = Random.default_rng()
    this_hop_rate, get_dt = hop_rates["simple_weak_strong_nr_hop_rate"]


    # setup parameters for model
    Ns = [num_particles_A, num_particles_B]
    num_spec = 2
    N = num_particles_A + num_particles_B
    L = system_size
    phys_t = phys_time
    dt = get_dt(num_particles_A, num_particles_B, alpha, chi)
    num_tsteps::Integer = floor(phys_time / dt)
    bc = bc
    alpha = alpha
    chi = chi
    hop_directions = (1, -1)

    # handle function input and create necessary parameters for everything
    model = ModelParameters(
        Ns, 
        num_spec,
        N,
        L,
        phys_t,
        num_tsteps,
        dt,
        bc,
        alpha,
        chi,
        hop_directions
    )
    state = NSpeciesState(model, rng)
    modelfunc = ModelFunctions(this_hop_rate, bcs[bc])
    chunk = Chunks(chunk_size, time_steps_between_snapshots, num_tsteps, dt)    

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
    return (model, state, chunk)

end
