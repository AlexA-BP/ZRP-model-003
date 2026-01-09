using HDF5
using ProgressBars

function mock_sim!(particles, prm)
    for i in 1:prm.N
        x = particles[i]
        particles[i] = mod(x-1 + rand([1, -1]), prm.L) + 1 
    end
end

function getlattice(particles, prm)
    lat = zeros(Int, prm.L)
    for x in particles
        lat[x] += 1
    end
    return lat
end


struct Parameters{T<:Integer}
    N::T
    L::T
    t::T
end


function main(N, L, t, chunkt)
    
    prm = Parameters(N, L, t)
    chunkN = prm.N
    chunkL = prm.L


    h5open("./data/test.h5", "w") do fid
        
        # setting up hdf5
        groupid = create_group(fid, "mocksim")
        dset_particle = create_dataset(fid, "mocksim/particles", datatype(Int), 
            dataspace((prm.N, prm.t)), chunk=(chunkN, chunkt))

        # dset_lattice = create_dataset(fid, "mocksim/lattice", datatype(Int), 
        #     dataspace((prm.t, prm.L)), chunk=(chunkt, chunkL))
        
        attrs(groupid)["N"] = prm.N
        attrs(groupid)["L"] = prm.L
        attrs(groupid)["t"] = prm.t


        ## setting up simulation
        particles = rand(1:prm.L, prm.N)
        # lat = getlattice(particles, prm)        

        ## setting up chunking
        # chunk_lat = Array{eltype(lat)}(undef, (prm.L, chunkt))

        chunk_particles = zeros(eltype(particles), (prm.N, chunkt))
        # chunk_lat = zeros(eltype(lat), (prm.L, chunkt))

        num_chunks = div(prm.t, chunkt)
        for ti in 1:num_chunks

            for tj in 1:chunkt
                chunk_particles[:, tj] = particles
                #chunk_lat[:, tj] = lat

                mock_sim!(particles, prm)
                # lat = getlattice(particles, prm)
            end

            dset_particle[:, ((ti-1)*chunkt + 1):(ti*chunkt)] = chunk_particles
            #dset_lattice[(ti-1)*chunkt + 1:ti*chunkt, :] = chunk_lat     
        end
    end
end
