using Random
const rng = Xoshiro(0)

using HDF5

include("sim/Structs.jl")
include("sim/BoundaryConditions.jl")
include("sim/kinetic_monte_carlo_single_species.jl")

using .BoundaryConditions

function setup_hdf5(
    fname::AbstractString, 
    params::Parameters,
    szrp::SingleZRP,
    chunk_size::Integer
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
            dataspace((params.N, params.t_tot));
            chunk=(params.N, chunk_size)
        )

        create_dataset(
            fid,
            "lattice",
            datatype(eltype(szrp.lattice)),
            dataspace((params.L, params.t_tot));
            chunk=(params.L, chunk_size)
        )
        
        create_dataset(
            fid,
            "times",
            datatype(Int),
            dataspace((params.t_tot,));
            chunk=(chunk_size,)
        )
    end
    return nothing
end

function run_and_write_chunked_simulation(
    fname::AbstractString,
    params::Parameters,
    szrp::SingleZRP,
    chunk_size::Integer
)
    h5open(fname, "r+") do fid 
        lattice_id = fid["lattice"]
        times_id = fid["times"]
        lattice_chunk = zeros(eltype(szrp.lattice), (params.L, chunk_size))
        times_chunk = zeros(Int, (chunk_size,))

        num_chunks = div(params.t_tot, chunk_size)
        for i in 1:num_chunks
            chunk_start = (i-1)*chunk_size
            chunk_end = i*chunk_size
            chunk_interval = (chunk_start+1):chunk_end

            fill_simulation_chunk!(
                lattice_chunk, 
                times_chunk, 
                params, 
                szrp, 
                chunk_size, 
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
    chunk_size::Integer,
    chunk_start::Integer,
)
    for ti in 1:chunk_size
        kinetic_monte_carlo_step!(
            szrp, params, chunk_start + ti
        )
        lattice_chunk[:, ti] = szrp.lattice
        times_chunk[ti] = chunk_start + ti
    end
    return nothing
end

function hop_rate(n::Integer, params::Parameters) 
    return 1. + params.b/n
end

function (@main)(
    fname::AbstractString,
    num_particles::Integer,
    system_size::Integer,
    bc::AbstractString,
    dt::Real,
    tot_timesteps::Integer,
    chunk_size::Integer,
    b::Real
) 

    # initialize simulation
    params = Parameters(
        num_particles,
        system_size,
        tot_timesteps,
        dt,
        bc,
        b,
    )
    szrp = SingleZRP(params, hop_rate, bcs[params.bc])


    # initialize HDF5 file
    setup_hdf5(fname, params, szrp, chunk_size,)

    # run simulation and write to HDF5 file
    run_and_write_chunked_simulation(fname, params, szrp, chunk_size,)

    return params, szrp
end
