using HDF5

struct Parameters{T<:Integer, S<:Real, U<:AbstractString}
    N::T
    L::T
    t_tot::T
    dt::S
    bc::U
end

struct SingleZRP{T<:Integer}
    particles::Vector{T}
    lattice::Vector{T}
end

function SingleZRP(params::Parameters)
    particles = _init_particles(params),
    lattice = _init_lattice(particles, params)
    return ZRP(particles, lattice)
end

function _init_particles(params::Parameters) 
    return rand(1:params.L, params.N)
end

function _init_lattice(particles, params::Parameters)
    lattice = zeros(eltype(particles), params.L)
    for x in particles
        lattice[x] += one(eltype(lattice))
    end
    return lattice
end


function setup_hdf5(params::Parameters, szrp::SingleZRP, chunksize::Integer)
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
            chunk=(params.N, chunksize)
        )

        create_dataset(
            fid,
            "lattice",
            datatype(eltype(szrp.lattice)),
            dataspace((params.L, params.t));
            chunk=(params.L, chunksize)
        )
    end
    return nothing
end


function (@main)(
    fname::AbstractString,
    num_particles::Integer,
    system_size::Integer,
    dt::Real,
    tot_timesteps::Integer,
    chunksize::Integer,
) 

    # initialize simulation
    prm = Parameters(
        num_particles,
        system_size,
        tot_timesteps,
        dt,
        bc,
    )
    szrp = SingleZRP(prm)

    # initialize HDF5 file
    setup_hdf5(params, szrp, chunksize)

    # run simulation and write to HDF5 file


end
