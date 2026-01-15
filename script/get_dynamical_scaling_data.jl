using HDF5
using ProgressBars

include("../src/julia/sim/utils.jl")
include("../src/julia/sim/basic_structs.jl")
include("../src/julia/sim/simulation.jl")

function (@main)(ARGS)

    N, L, t, chunk_t, id = ARGS
    N, L, t, chunk_t = parse.(Int, [N, L, t, chunk_t])

    #TODO make this actually good and not terrible

    dt = 1/4.5
    bc = "p"
    num_species = 1

    target_dir = "/data/lmp/aarnhold/dynamical_scaling/"
    target = "L_$(L)_$(id).h5"

    prm = Parameters(N, L, t, num_species, dt, bc)
    zrp = ZRP(prm)
    chunk = Chunk(zrp, prm, chunk_t)

    h5open(target_dir * target, "w") do fid
        
        group_id = create_group(fid, "sim")
        attrs(group_id)["N"] = prm.N
        attrs(group_id)["L"] = prm.L
        attrs(group_id)["t"] = prm.t
        attrs(group_id)["bc"] = prm.bc
        
        lattice_id = create_dataset(
            fid,
            "sim/lattice",
            datatype(eltype(zrp.lattice)),
            dataspace((prm.L*prm.num_species, prm.t)),
            chunk=(prm.L*prm.num_species, chunk.t)
        )
        
        for ti in 1:chunk.num
            for tj in 1:chunk.t
                assign_hyperslab!(chunk.particles, zrp.particles, tj)
                assign_hyperslab!(chunk.lattice, zrp.lattice, tj)

                update!(zrp, prm, periodic_bc)
            end
            ts = ((ti-1)*chunk.t + 1):(ti*chunk.t)
            lattice_id[:, ts] = chunk.lattice
        end
    end

    println("Finished L $(L) and $(id)")

    return nothing
end