using HDF5

function maxs(h5filepath, t_skip)
    h5open(h5filepath, "r+") do fid
        dset = read_dataset(fid, "group_test/lattice")               
        t = read_attribute(fid["group_test"], "t")        
        maxs = Vector{Int}(undef, div(t, t_skip))
        for (i, col) in enumerate(eachcol(dset[begin:t_skip:end]))
            maxs[i] = maximum(col)
        end
        
        group_id = create_group(fid, "post_processing")
        group_id["maximums"] = maxs
        return nothing
    end
end

function prob_site_has_n_particles(h5filepath)
    h5open(h5filepath, "r+") do fid
        dset = read_dataset(fid, "group_test/lattice")
        N = read_attribute(fid["group_test"], "N")
        t = read_attribute(fid["group_test"], "t")

        prob_site_n = zeros(N)

        for lattice in eachcol(dset)
            for occupation in eachindex(lattice)
                prob_site_n[occupation] += 1.
            end
        end

        prob_site_n /= t

        fid["post_processing/prob_site_has_n_particles"] = prob_site_n
        return nothing

    end
end
