using HDF5
using ProgressBars
using Distributions

include("utils.jl")
include("basic_structs.jl")
include("simulation.jl")

function main(N, L, t, num_species, dt, bc, chunk_t, fname="./data/test.h5")

    us = ones(L)
    
    us[rand(1:L, 2)] .= 0.2 
    dt = 1/minimum(us)

    prm = Parameters(N, L, t, num_species, dt, bc)
    zrp = ZRP(prm)
    chunk = Chunk(zrp, prm, chunk_t)

    h5open(fname, "w") do fid
        
        attrs(fid)["N"] = prm.N
        attrs(fid)["L"] = prm.L
        attrs(fid)["t"] = prm.t
        attrs(fid)["bc"] = prm.bc
        
        particles_id = create_dataset(
            fid,
            "particles",
            datatype(eltype(zrp.particles)),
            dataspace((prm.N, prm.t)),
            chunk=(prm.N, chunk.t)
        )

        lattice_id = create_dataset(
            fid,
            "lattice",
            datatype(eltype(zrp.lattice)),
            dataspace((prm.L*prm.num_species, prm.t)),
            chunk=(prm.L*prm.num_species, chunk.t)
        )
        
        hoprates_id = create_dataset(
            fid,
            "hoprates",
            datatype(Float64),
            dataspace((prm.L*prm.num_species, 1)),
        )

        hoprates_id[:, 1] = us

        for ti in ProgressBar(1:chunk.num)
            for tj in 1:chunk.t
                assign_hyperslab!(chunk.particles, zrp.particles, tj)
                assign_hyperslab!(chunk.lattice, zrp.lattice, tj)

                update!(zrp, prm, periodic_bc, us)
            end
            ts = ((ti-1)*chunk.t + 1):(ti*chunk.t)
            particles_id[:, ts] = chunk.particles
            lattice_id[:, ts] = chunk.lattice
        end
    end
    return nothing
end