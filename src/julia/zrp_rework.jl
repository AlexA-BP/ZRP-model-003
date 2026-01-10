using HDF5

include("sim/utils.jl")
include("sim/basic_structs.jl")
include("sim/simulation.jl")

function main(N, L, t, num_species, dt, bc, chunk_t)

    prm = Parameters(N, L, t, num_species, dt, bc)
    zrp = ZRP(prm)
    chunk = Chunk(zrp, prm, chunk_t)

    # return prm, zrp, chunk
    h5open("./data/test.h5", "w") do fid
        
        group_id = create_group(fid, "group_test")
        attrs(group_id)["N"] = prm.N
        attrs(group_id)["L"] = prm.L
        attrs(group_id)["t"] = prm.t
        
        particles_id = create_dataset(
            fid,
            "group_test/particles",
            datatype(eltype(zrp.particles)),
            dataspace((prm.N, prm.t)),
            chunk=(prm.N, chunk.t)
        )

        lattice_id = create_dataset(
            fid,
            "group_test/lattice",
            datatype(eltype(zrp.lattice)),
            dataspace((prm.L, prm.num_species, prm.t)),
            chunk=(prm.L, prm.num_species, chunk.t)
        )
        
        for ti in 1:chunk.num
            for tj in 1:chunk.t
                assign_hyperslab!(chunk.particles, zrp.particles, tj)
                assign_hyperslab!(chunk.lattice, zrp.lattice, tj)

                # update!(zrp, prm)
            end
            # indeces = ((ti-1)*chunk.t + 1):(ti*chunk.t)
            # particles_id[:, indeces] = chunk.particles
            # lattice_id[:, :, indeces] = chunk.lattice
        end
    end
    return nothing
end