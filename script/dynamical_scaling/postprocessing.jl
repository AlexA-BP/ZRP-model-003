using HDF5
using Statistics
using ProgressBars

const src_dir = "/data/lmp/aarnhold/dynamical_scaling/"

function get_ensemble_maxs(sim_name, chunksize) 
    
    # get files
    files = readdir(joinpath(src_dir,sim_name))
    sims = joinpath.(src_dir, sim_name, files)
    nensemble = length(sims)

    # setting up parameters for later use
    keys = ["L", "N", "t", "bc"]
    parameters = get_parameters(sims[begin], keys)

    # get ensemble maxs
    nchunks = div(parameters["t"], chunksize)
    h5open(joinpath(src_dir, sim_name*".h5"), "cw") do fid

        write_parameters(fid, parameters)

        ensemble_maxs_id = create_dataset(
            fid, 
            "ensemble_maxs",
            datatype(Int),
            dataspace(parameters["t"], nensemble),
            chunk=(chunksize, nensemble)
        )
        calculate_maxs_time_series_ensemble!(
            ensemble_maxs_id, 
            sims, 
            chunksize, 
            nchunks
        )
    end
    return nothing
end

function get_mean_maxs(sim_name, chunksize) 
    
    
    h5open(joinpath(src_dir, sim_name*".h5"), "cw") do fid

        nchunks = div(read_attribute(fid, "t"), chunksize)
        ensemble_maxs_id = open_dataset(fid, "ensemble_maxs")

        dset_name = "avg_maxs"
        mean_maxs_id = create_dataset(
            fid, 
            dset_name,
            datatype(Float64),
            dataspace(read_attribute(fid, "t"), 1),
            chunk=(chunksize, 1)
        )

        calculate_maxs_time_series_mean!(
            mean_maxs_id,
            ensemble_maxs_id,
            chunksize,
            nchunks
        )

    end
    return nothing
end

function calculate_maxs_time_series_ensemble!(
    ensemble_maxs_id, 
    sims, 
    chunksize, 
    nchunks
)
    for (i, sim) in ProgressBar(enumerate(sims))
        h5open(sim, "r") do fid
            calculate_maxs_time_series!(ensemble_maxs_id, i, fid, chunksize, nchunks)
        end
    end
    return nothing
end

function calculate_maxs_time_series!(
    ensemble_maxs_id, 
    ensemble_index, 
    subfid, 
    chunksize, 
    nchunks
)
    maxs = zeros(Int, chunksize)
    for n in 1:nchunks
        ts = ((n-1)*chunksize + 1):(n*chunksize)
        ensemble_lattices = open_dataset(subfid, "sim/lattice")
        ensemble_lattices_chunk = ensemble_lattices[:, ts]
        for (t, lattice) in enumerate(eachcol(ensemble_lattices_chunk))
            maxs[t] = maximum(lattice)
        end
        ensemble_maxs_id[ts, ensemble_index] = maxs                    
    end
    return nothing
end

function calculate_maxs_time_series_mean!(
    mean_maxs_id,
    ensemble_maxs_id, 
    chunksize,
    nchunks)
    for n in 1:nchunks
        ts = ((n-1)*chunksize+1):(n*chunksize)
        mean_maxs_id[ts, 1] = mean(ensemble_maxs_id[ts, :], dims=2)
    end
    return nothing
end

function get_parameters(sim, keys)
    parameters = Dict()
    h5open(sim, "r") do fid
        for key in keys
            parameters[key] = read_attribute(fid["sim"], key)
        end
    end
    return parameters
end

function write_parameters(fid, parameters)
    for (key, value) in parameters
        attributes(fid)[key] = value
    end
end