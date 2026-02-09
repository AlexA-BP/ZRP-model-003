using Random
const rng = Xoshiro(0)

include("sim/Structs.jl")
include("sim/BoundaryConditions.jl")
include("sim/kinetic_monte_carlo_single_species.jl")
include("sim/filehandling.jl")

using .BoundaryConditions

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
        num_particles, system_size, tot_timesteps, dt, bc, b,
    )
    szrp = SingleZRP(params, hop_rate, bcs[params.bc])

    # initialize HDF5 file
    setup_hdf5(fname, params, szrp, chunk_size,)

    # run simulation and write to HDF5 file
    run_and_write_chunked_simulation(fname, params, szrp, chunk_size,)

    return nothing
end
