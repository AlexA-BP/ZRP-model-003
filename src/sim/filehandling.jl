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

        attrs(fid)["N_A"] = model.Ns[1]
        attrs(fid)["N_B"] = model.Ns[2]
        attrs(fid)["N"] = model.N
        attrs(fid)["L"] = model.L
        attrs(fid)["rho_A"] = model.Ns[1] / model.L
        attrs(fid)["rho_B"] = model.Ns[2] / model.L
        attrs(fid)["rho"] = model.N / model.L
        attrs(fid)["t_tot"] = model.t_tot
        attrs(fid)["dt"] = model.dt
        attrs(fid)["bc"] = model.bc
        attrs(fid)["alpha"] = model.alpha
        attrs(fid)["chi"] = model.chi
        attrs(fid)["chunk_size"] = chunk.size
        attrs(fid)["chunk_num"] = chunk.num
        attrs(fid)["chunk_time_steps_between_snapshots"] = chunk.time_steps_between_saves


        fid["species"] = state.species
        
        create_dataset(
            fid,
            "particles",
            datatype(eltype(state.particles)),
            dataspace((model.N, chunk.num*chunk.size)),
            chunk=(model.N, chunk.size),
            blosc=3,
        )

        create_dataset(
            fid,
            "lattice",
            datatype(eltype(state.lattice)),
            dataspace((model.L, model.num_spec, chunk.num*chunk.size)),
            chunk=(model.L, model.num_spec, chunk.size),
            blosc=3,
        )
        

        create_dataset(
            fid,
            "times",
            datatype(Int),
            dataspace((chunk.num*chunk.size,)),
            chunk=(chunk.size,),
            blosc=3,
        )
    end
    return nothing
end

function run_and_write_chunked_simulation(
    fname::AbstractString,
    model::ModelParameters,
    state::NSpeciesState,
    modelfunc::ModelFunctions,
    chunk::Chunks,
    rng::AbstractRNG,
)
    h5open(fname, "r+") do fid 
        lattice_id = fid["lattice"]
        times_id = fid["times"]

        lattice_chunk = zeros(
            eltype(state.lattice), (model.L, model.num_spec, chunk.size)
        )
        times_chunk = zeros(Int, (chunk.size,))

        for i in 1:chunk.num
            chunk_start = (i-1)*chunk.size
            chunk_end = i*chunk.size
            chunk_interval = (chunk_start+1):chunk_end

            fill_simulation_chunk!(
                lattice_chunk, 
                times_chunk, 
                model, 
                state, 
                modelfunc,
                chunk,
                chunk_start,
                rng,
            )
    
            # lattice_id[:, :, chunk_interval] = lattice_chunk
            # times_id[chunk_interval] = times_chunk
        end
    end
    return nothing
end

function fill_simulation_chunk!(
    lattice_chunk::AbstractArray,
    times_chunk::AbstractVector,
    model::ModelParameters,
    state::NSpeciesState,
    modelfunc::ModelFunctions,
    chunk::Chunks,
    chunk_start::Integer,
    rng::AbstractRNG,
)
    for k in 1:chunk.size
        for ti in 1:chunk.time_steps_between_saves
            kinetic_monte_carlo_step!(
                state, 
                model, 
                modelfunc,
                (chunk_start + (k-1))*chunk.time_steps_between_saves + ti,
                rng,
            )
        end
        # lattice_chunk[:, :, k] = state.lattice
        # times_chunk[k] = (chunk_start + k)*chunk.time_steps_between_saves
    end
    return nothing
end