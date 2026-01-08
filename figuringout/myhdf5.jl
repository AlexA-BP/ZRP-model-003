using Revise
using HDF5

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


struct Parameters
    N::Integer
    L::Integer
    t::Integer
end


function main(N, L, t, chunkt)
    
    prm = Parameters(N, L, t)

    h5open("./data/test.h5", "w") do fid
        
        # setting up hdf5
        groupid = create_group(fid, "mocksim")
        dset_particle = create_dataset(fid, "mocksim/particles", datatype(Int), 
            dataspace((prm.t, prm.N)), chunk=(chunkt, chunkN))

        dset_lattice = create_dataset(fid, "mocksim/lattice", datatype(Int), 
            dataspace((prm.t, prm.L)), chunk=(chunkt, chunkL))
        
        attrs(groupid)["N"] = prm.N
        attrs(groupid)["L"] = prm.L
        attrs(groupid)["t"] = prm.t

        # setting up chunking
        chunk_particles = Array{eltype(particles)}(undef, (chunkt, prm.N))
        chunk_lattice = Array{eltype(lattice)}(undef, (chunkt, prm.L))
                
        # setting up simulation
        particles = [rand(1:prm.L) for _ = 1:prm.N]
        lat = getlattice(particles, prm)        

        for ti in 1:div(prm.t, chunkt)

            for tj in 1:chunkt
                chunk_particles[:, tj] = particles
                chunk_lattice[:, tj] = lattice
                
                mock_sim!(particles, prm)
                lat = getlattice(particles, prm)
            end

            # TODO: get the indeces right
            dset_particle[ti*chunkt:(ti+1)*chunkt, :] = chunk_particles
            dset_lattice[ti, :] = lat      
        end
        

    end

end
