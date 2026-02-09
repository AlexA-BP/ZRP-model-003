using HDF5

function setup_hdf5(
    fname::AbstractString, 
    params::Parameters,
    szrp::SingleZRP,
    chunk::Chunks,
)
    h5open(fname, "w") do fid

        attrs(fid)["N"] = params.N
        attrs(fid)["L"] = params.L
        attrs(fid)["t_tot"] = params.t_tot
        attrs(fid)["dt"] = params.dt
        attrs(fid)["bc"] = params.bc

        create_dataset(
            fid,
            "particles",
            datatype(eltype(szrp.particles)),
            dataspace((params.N, chunk.num*chunk.size));
            chunk=(params.N, chunk.size)
        )

        create_dataset(
            fid,
            "lattice",
            datatype(eltype(szrp.lattice)),
            dataspace((params.L, chunk.num*chunk.size));
            chunk=(params.L, chunk.size)
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
    params::Parameters,
    szrp::SingleZRP,
    chunk::Chunks,
)
    h5open(fname, "r+") do fid 
        lattice_id = fid["lattice"]
        times_id = fid["times"]
        lattice_chunk = zeros(eltype(szrp.lattice), (params.L, chunk.size))
        times_chunk = zeros(Int, (chunk.size,))

        for i in 1:chunk.num
            chunk_start = (i-1)*chunk.size
            chunk_end = i*chunk.size
            chunk_interval = (chunk_start+1):chunk_end

            fill_simulation_chunk!(
                lattice_chunk, 
                times_chunk, 
                params, 
                szrp, 
                chunk,
                chunk_start,
            )
    
            lattice_id[:, chunk_interval] = lattice_chunk
            times_id[chunk_interval] = times_chunk
        end
    end
    return nothing
end

function fill_simulation_chunk!(
    lattice_chunk::AbstractArray,
    times_chunk::AbstractVector,
    params::Parameters,
    szrp::SingleZRP,
    chunk::Chunks,
    chunk_start::Integer,
)
    for k in 1:chunk.size
        for ti in 1:chunk.saving_time_step
            kinetic_monte_carlo_step!(
                szrp, 
                params, 
                chunk_start + (k-1)*chunk.saving_time_step + ti
            )
        end
        lattice_chunk[:, k] = szrp.lattice
        times_chunk[k] = chunk_start + k*chunk.saving_time_step
    end
    return nothing
end