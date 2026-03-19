#=
Czech-list
    [ ] basic structure
    [ ] file handling
        [ ]
        [ ]
    [ ] simulation  
=#
using Random
using Format

include("./sim/TwoSpeciesZRP/Structs.jl")
include("./sim/TransitionRates.jl")
include("./sim/BoundaryConditions.jl")
include("./sim/filehandling.jl")
include("./sim/NaiveMonteCarlo.jl")
include("./sim/RejectionKineticMonteCarlo.jl")

import .TransitionRates: TransitionRatesDict
using .BoundaryConditions


function main_dev(;
    # function paramters, for now hardcoded
    num_particles_A=5,
    num_particles_B=5,
    system_size=10,
    bc="p",
    alpha=0., 
    chi=3., 
    phys_time=10.,
    chunk_size=10,
    chunk_num=10,
    snapshot_phys_time_diff=1,
    dirname="./data/twospecies/"
)
    

    rng = Xoshiro(2) 
    # rng = Random.default_rng()
    this_hop_rate = (
        TransitionRatesDict["simple_weak_strong_nr_hop_rate"]
    )


    # setup parameters for model
    Ns = [num_particles_A, num_particles_B]
    num_spec = 2
    N = num_particles_A + num_particles_B
    L = system_size
    phys_t = phys_time
    bc = bc
    alpha = alpha
    chi = chi
    hop_directions = (1, -1)


    # handle function input and create necessary parameters for everything
    model = ModelParameters(
        Ns, 
        L,
        phys_t,
        num_spec,
        N,
        bc,
        alpha,
        chi,
        hop_directions
    )

    modelfunc = ModelFunctions(this_hop_rate, bcs[bc])
    state = NSpeciesState(model, modelfunc, rng)
    # chunk = NaiveChunks(chunk_size, time_steps_between_snapshots, num_tsteps, dt)    
    chunk = Chunks(chunk_size, chunk_num, snapshot_phys_time_diff)    


    # set up file structure
    if !isdir(dirname)
        mkdir(dirname) 
    end
    fcnt = 1
    while state.t <= model.phys_t
        
        fname = format(dirname * "trial_{1:03d}", fcnt)
        # setup hdf5
        setup_hdf5(fname, model, state, chunk)
        
        # return (fname, model, state, modelfunc, chunk, rng)
        run_and_write_chunked_simulation(
            fname,
            model,
            state,
            modelfunc,
            chunk,
            rng,
        )

        fcnt += 1
    end
    # return (model, state, chunk)
    return nothing


end
