using Revise
using Printf


struct Parameters
    L::Integer
    N::Integer
    b::Real
    dt::Real
    bc::String
end

function get_lattice(particles, prm::Parameters)
    lat = zeros(prm.L)
    for pos in particles
        lat[pos] += 1
    end
    return lat
end

function update_particle_and_lattice!(particles, lat, i, x, step, prm::Parameters)
    new_x = apply_bc!(x + step, prm)
    particles[i] = new_x
    lat[x] -= 1
    lat[new_x] += 1
    return nothing
end

function update_particle!(particles, i, x, step, prm::Parameters)1
    new_x = apply_bc!(x + step, prm)
    particles[i] = new_x
end

function update_lattice!(lat, x, step, prm::Parameters)
    new_x = apply_bc!(x + step, prm)
    lat[x] -= 1
    lat[new_x] += 1
    return nothing
end

function apply_bc!(x, prm::Parameters)
    if prm.bc in ["p", "periodic"]
        return apply_bc_periodic!(x, prm)
    else
        throw("$(prm.bc) is a not a supported bc!")
    end
    return nothing
end
        

function apply_bc_periodic!(x, prm::Parameters)
    if x > prm.L
        x -= prm.L
    elseif x <= 0
        x += prm.L
    end
    return x
end

function select_sim(updating_scheme)
    if updating_scheme == ("particle", "sequential")
        return sim_particlebased_sequential
    elseif updating_scheme == ("particle", "parallel")
        return sim_particlebased_parallel
    elseif updating_scheme == ("lattice", "sequential")
        return sim_latticebased_sequential
    elseif updating_scheme == ("lattice", "parallel")
        return sim_latticebased_parallel
    end
end

function sim_particlebased_parallel(particles, lat, prm::Parameters)
    for i = 1:prm.N
        x = particles[i]
        n = lat[x]

        prob = hop_rate(n, prm)*prm.dt/n

        if rand() < prob
            step = get_step(x)
            update_particle!(particles, i, x, step, prm)
        end
    end
    lat = get_lattice(particles, prm)
    return particles, lat
end

function sim_particlebased_sequential(particles, lat, prm)
    for _ = 1:prm.N
        randi = rand(1:prm.N)
        x = particles[randi]
        n = lat[x]

        prob = hop_rate(n, prm)*prm.dt/n

        if rand() < prob
            step = get_step(x)
            update_particle_and_lattice!(particles, lat, randi, x, step, prm)
        end
    end
    return particles, lat
end
        
function sim_latticebased_sequential(particles, lat, prm::Parameters)
    for _ = 1:prm.L
        randx = rand(1:prm.L)
        n = lat[randx]
        if n != 0
            prob = hop_rate(n, prm)*prm.dt
            if rand() < prob
                step = get_step(randx)
                update_lattice!(lat, randx, step, prm)
            end
        end
    end
    return particles, lat
end

function sim_latticebased_parallel(particles, lat, prm::Parameters)
    lat_temp = copy(lat)
    for x in 1:prm.L
        n = lat_temp[x]
        if n != 0
            prob = hop_rate(n, prm)*prm.dt
            if rand() < prob
                step = get_step(x)
                update_lattice!(lat, x, step, prm)
            end
        end
    end
    return particles, lat
end

function get_step(x)
    return oneunit(x)
end

function hop_rate(n, prm::Parameters)
    return (1 + prm.b/n)
end

function main(;L, N, b, t, bc="p", updating_scheme=("particle", "sequential"))

    dt = 1/(1 + b)
    t = trunc(t/dt)      

    prm = Parameters(L, N, b, dt, bc)       

    particles = [rand(1:L) for _ in 1:N]
    lat = zeros(L)
    for pos in particles
        lat[pos] += 1
    end

    sim = select_sim(updating_scheme)

    open("test.dat", "w") do f
        for i in 1:t

            if mod(i,100) == 0
                println(f, lat)
            end
            
            particles, lat = sim(particles, lat, prm)
            
            @printf("progress %.2f%%\r", 100*i/t)

        end
    end

    return lat

end
