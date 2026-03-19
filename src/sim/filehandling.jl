using HDF5
using H5Zblosc
using ProgressBars

"""
Set up hdf5 file. Save parameters and chunks as attributes
and create the datasets.   
"""
function setup_hdf5(
    fname::AbstractString, 
    model::ModelParameters,
    state::NSpeciesState,
    chunk::Chunks,
)
    h5open(fname, "w") do fid


        # save all parameters to hdf5 as attributes
        attrs(fid)["N_A"] = model.Ns[1]
        attrs(fid)["N_B"] = model.Ns[2]
        attrs(fid)["N"] = model.N
        attrs(fid)["L"] = model.L
        attrs(fid)["rho_A"] = model.Ns[1] / model.L
        attrs(fid)["rho_B"] = model.Ns[2] / model.L
        attrs(fid)["rho"] = model.N / model.L
        attrs(fid)["phys_time"] = model.phys_t
        attrs(fid)["bc"] = model.bc
        attrs(fid)["alpha"] = model.alpha
        attrs(fid)["chi"] = model.chi
        attrs(fid)["chunk_size"] = chunk.size
        attrs(fid)["chunk_num"] = chunk.num
        attrs(fid)["chunk_snapshot_phys_time_diff"] = (
            chunk.snapshot_phys_time_diff)


        # initialize the datasets
        blosc_compression_level = 3


        # species can be saved immediately, since it stays constant throughout
        # the simulation
        fid[
            "species", 
            chunk=(model.N,), 
            blosc=blosc_compression_level
        ] = state.species
        
        # initialize empty "particles", "lattice" and "times" datasets

        create_dataset(
            fid,
            "particles",
            datatype(eltype(state.particles)),
            dataspace((model.N, chunk.num*chunk.size)),
            chunk=(model.N, chunk.size),
            blosc=blosc_compression_level ,
        )

        create_dataset(
            fid,
            "lattice",
            datatype(eltype(state.lattice)),
            dataspace((model.L, model.num_spec, chunk.num*chunk.size)),
            chunk=(model.L, model.num_spec, chunk.size),
            blosc=blosc_compression_level,
        )
        

        create_dataset(
            fid,
            "times",
            datatype(Float64),
            dataspace((chunk.num*chunk.size,)),
            chunk=(chunk.size,),
            blosc=blosc_compression_level,
        )
    end
    return nothing
end

"""
Run the simulation in chunks and write the chunks to HDF5. 
"""
function run_and_write_chunked_simulation(
    fname::AbstractString,
    model::ModelParameters,
    state::NSpeciesState,
    modelfunc::ModelFunctions,
    chunk::Chunks,
    rng::AbstractRNG,
)
    h5open(fname, "r+") do fid 

        # open datasets for saving
        particles_id = fid["particles"]
        lattice_id = fid["lattice"]
        times_id = fid["times"]

        # init in-memory chunks for saving data
        particles_chunk = zeros(eltype(state.particles), (model.N, chunk.size))
        lattice_chunk = zeros(
            eltype(state.lattice), (model.L, model.num_spec, chunk.size)
        )
        times_chunk = zeros(Float64, (chunk.size,))

        # loop over chunks
        for i in 1:chunk.num

            # 
            chunk_start_time = (i-1)*chunk.size
            chunk_end_time = i*chunk.size
            chunk_time_interval = (chunk_start_time+1):chunk_end_time

            fill_simulation_chunk_rejection_kinetic_monte_carlo!(
                particles_chunk, 
                lattice_chunk, 
                times_chunk, 
                model, 
                state, 
                modelfunc,
                chunk,
                chunk_start_time,
                rng,
            )
    
            particles_id[:, chunk_time_interval] = particles_chunk
            lattice_id[:, :, chunk_time_interval] = lattice_chunk
            times_id[chunk_time_interval] = times_chunk
        end
    end
    return nothing
end

"""
Fill the simulation chunks for "Particles", "Lattice" and "Times" for the 
naive monte carlo algorithm.  
"""
function fill_simulation_chunk_naive_monte_carlo!(
    particles_chunk::AbstractArray,
    lattice_chunk::AbstractArray,
    times_chunk::AbstractVector,
    model::ModelParameters,
    state::NSpeciesState,
    modelfunc::ModelFunctions,
    chunk::Chunks,
    chunk_start_time::Integer,
    rng::AbstractRNG,
)
    for k in 1:chunk.size
        for ti in 1:chunk.snapshot_phys_time_diff

            current_time = (
                (chunk_start_time + (k-1))*chunk.snapshot_phys_time_diff + ti
            )

            naive_monte_carlo_step!(
                state, 
                model, 
                modelfunc,
                current_time,
                rng,
            )
        end

        particles_chunk[:, k] = state.particles
        lattice_chunk[:, :, k] = state.lattice
        times_chunk[k] = (chunk_start_time + k)*chunk.snapshot_phys_time_diff
    end
    return nothing
end

"""
Fill the simulation chunks for "Particles", "Lattice" and "Times" for the 
rejection kinetic monte carlo algorithm.  
"""
function fill_simulation_chunk_rejection_kinetic_monte_carlo!(
    particles_chunk::AbstractArray,
    lattice_chunk::AbstractArray,
    times_chunk::AbstractVector,
    model::ModelParameters,
    state::NSpeciesState,
    modelfunc::ModelFunctions,
    chunk::Chunks,
    chunk_start_time::Integer,
    rng::AbstractRNG,
)
    for k in 1:chunk.size
        current_phys_time_diff = 0
        while current_phys_time_diff < chunk.snapshot_phys_time_diff

            rejection_kinetic_monte_carlo_step!(
                state, 
                model, 
                modelfunc,
                rng,
            )

            current_phys_time_diff += state.dt
            
        end


        particles_chunk[:, k] = state.particles
        lattice_chunk[:, :, k] = state.lattice
        times_chunk[k] = state.t
    end
    return nothing
end