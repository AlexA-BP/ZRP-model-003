using Revise
using ProgressBars
using Statistics
using GLMakie

includet("./sim/Simulations.jl")
includet("./condensatescaling/CondensateScaling.jl")
using .Simulations
using .CondensateScaling

function simulate(;L, N, b, t, bc="p", updating_scheme=("particle", "sequential"))
    dt = 1/(1 + b)
    t = trunc(t/dt)      

    prm = Parameters(L, N, b, dt, bc)       

    particles = init_particles(prm)
    lat = get_lattice(particles, prm)
    sim = selectsim(updating_scheme)

    open("data/test.dat", "w") do f
        for i in ProgressBar(1:t)

            if mod(i-1,100) == 0
                println(f, lat)
            end
            particles, lat = sim(particles, lat, prm)

        end
    end
    return lat
end

function condensatescaling(
    ;
    L::Integer, 
    N::Integer, 
    b::Real, 
    t::Real, 
    bc::String="p", 
    updating_scheme::Tuple{String, String}=("particle", "sequential"),
    nsamples::Integer,
)
    dt = 1/(1 + b)
    tsteps::Integer = trunc(t/dt)      

    prm = Parameters(L, N, b, dt, bc)       

    sim! = selectsim(updating_scheme)
    particles = init_particles(prm)
    lat = get_lattice(particles, prm)

    mean_condensate_sizes = zeros(tsteps, nsamples)
    for sample in ProgressBar(1:nsamples)
        for tstep in 1:tsteps
            mean_condensate_sizes[tstep, sample] =  mean_condensate_size(lat, 3)
            particles, lat = sim!(particles, lat, prm)
        end
    end

    final_result = mean(mean_condensate_sizes, dims=2)

    # fig, ax, p = lines(collect(1:tsteps), final_result)

    return mean_condensate_sizes
end