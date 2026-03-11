#=
Czech-list
    [ ] basic structure
    [ ] file handling
        [ ]
        [ ]
    [ ] simulation  
=#
using Random
const rng = Xoshiro(2) 



function main_dev()
    # function paramters, for now hardcoded
    num_particles_A = 5
    num_particles_B = 6
    system_size = 10
    bc = "p"
    tot_timesteps = 100
    chunk_size = 10
    time_steps_between_snapshots = 5
    
    # handle function input and create necessary parameters for everything
    prm = TwoSpeciesParameters()
    sim = TwoSpeciesZRP()
    chunk = Chunks(chunk_size, time_steps_between_snapshots, tot_timesteps)    

    # setup hdf5


    # perform simulation and write data


end