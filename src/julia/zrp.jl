using Revise
using ProgressBars

includet("./sim/Simulations.jl")
using .Simulations

function (@main)(;L, N, b, t, bc="p", updating_scheme=("particle", "sequential"))

    dt = 1/(1 + b)
    t = trunc(t/dt)      

    prm = Parameters(L, N, b, dt, bc)       

    particles = init_particles(prm)
    lat = get_lattice(particles, prm)

    sim = selectsim(updating_scheme)

    open("data/test.dat", "w") do f
        for i in ProgressBar(1:t)

            if mod(i,100) == 0
                println(f, lat)
            end
            
            particles, lat = sim(particles, lat, prm)
            
            #@printf("progress %.2f%%\r", 100*i/t)

        end
    end

    return lat

end
