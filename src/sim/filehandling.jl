using HDF5
using ProgressBars

"""
Set up hdf5 file. 
"""
function setup_hdf5(
    fname::AbstractString, 
    prm::TwoSpeciesParameters,
    state::TwoSpeciesState,
    chunk::Chunks,
)
    h5open(fname, "w") do fid

        attrs(fid)["N_A"] = prm.N_A
        attrs(fid)["N_B"] = prm.N_B
        attrs(fid)["N"] = prm.N
        attrs(fid)["L"] = prm.L
        attrs(fid)["rho_A"] = prm.N_A / prm.L
        attrs(fid)["rho_B"] = prm.N_B / prm.L
        attrs(fid)["rho"] = prm.N / prm.L
        attrs(fid)["t_tot"] = prm.t_tot
        attrs(fid)["dt"] = prm.dt
        attrs(fid)["bc"] = prm.bc
        attrs(fid)["alpha"] = prm.alpha
        attrs(fid)["chi"] = prm.chi
        attrs(fid)["chunk_size"] = chunk.size
        attrs(fid)["chunk_num"] = chunk.num
        attrs(fid)["chunk_time_steps_between_snapshots"] = chunk.time_steps_between_saves

        create_dataset(
            fid,
            "particles_A",
            datatype(eltype(state.particles_A)),
            dataspace((prm.N_A, chunk.num*chunk.size));
            chunk=(prm.N_A, chunk.size)
        )

        create_dataset(
            fid,
            "lattice_A",
            datatype(eltype(state.lattice_A)),
            dataspace((prm.L, chunk.num*chunk.size));
            chunk=(prm.L, chunk.size)
        )
        
        create_dataset(
            fid,
            "particles_B",
            datatype(eltype(state.particles_B)),
            dataspace((prm.N_B, chunk.num*chunk.size));
            chunk=(prm.N_B, chunk.size)
        )

        create_dataset(
            fid,
            "lattice_B",
            datatype(eltype(state.lattice_B)),
            dataspace((prm.L, chunk.num*chunk.size));
            chunk=(prm.L, chunk.size)
        )

        create_dataset(
            fid,
            "times",
            datatype(Int),
            dataspace((chunk.num*chunk.size,));
            chunk=(chunk.size,)
        )
    end
    return nothing
end

function run_and_write_chunked_simulation(
    fname::AbstractString,
    prm::Parameters,
    state::TwoSpeciesState,
    modelfunc::ModelFunctions,
    chunk::Chunks,
)
    h5open(fname, "r+") do fid 
        lattice_A_id = fid["lattice_A"]
        lattice_B_id = fid["lattice_B"]
        times_id = fid["times"]

        lattice_A_chunk = zeros(eltype(state.lattice_A), (prm.L, chunk.size))
        lattice_B_chunk = zeros(eltype(state.lattice_B), (prm.L, chunk.size))
        times_chunk = zeros(Int, (chunk.size,))

        for i in 1:chunk.num
            chunk_start = (i-1)*chunk.size
            chunk_end = i*chunk.size
            chunk_interval = (chunk_start+1):chunk_end

            fill_simulation_chunk!(
                lattice_A_chunk, 
                lattice_B_chunk, 
                times_chunk, 
                prm, 
                state, 
                modelfunc,
                chunk,
                chunk_start,
            )
    
            lattice_A_id[:, chunk_interval] = lattice_A_chunk
            lattice_B_id[:, chunk_interval] = lattice_B_chunk
            times_id[chunk_interval] = times_chunk
        end
    end
    return nothing
end

function fill_simulation_chunk!(
    lattice_A_chunk::AbstractArray,
    lattice_B_chunk::AbstractArray,
    times_chunk::AbstractVector,
    prm::Parameters,
    state::TwoSpeciesState,
    modelfunc::ModelFunctions,
    chunk::Chunks,
    chunk_start::Integer,
)
    for k in 1:chunk.size
        for ti in 1:chunk.time_steps_between_saves
            kinetic_monte_carlo_step!(
                state, 
                prm, 
                modelfunc,
                (chunk_start + (k-1))*chunk.time_steps_between_saves + ti
            )
        end
        lattice_A_chunk[:, k] = state.lattice_A
        lattice_B_chunk[:, k] = state.lattice_B
        times_chunk[k] = (chunk_start + k)*chunk.time_steps_between_saves
    end
    return nothing
end