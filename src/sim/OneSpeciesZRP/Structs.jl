struct Parameters{T<:Integer, S<:Real, U<:AbstractString}
    N::T
    L::T
    t_tot::T
    dt::S
    bc::U
    b::S
end

abstract type ZRP end

struct SingleZRP{T<:Integer, F<:Function, G<:Function}
    particles::Vector{T}
    lattice::Vector{T}
    hop_rate::F
    bc::G
end

function SingleZRP(
    params::Parameters, 
    hop_rate::F, 
    bc::G
) where {F<:Function, G<:Function}
    particles = _init_particles(params)
    lattice = _init_lattice(particles, params)
    return SingleZRP(particles, lattice, hop_rate, bc)
end

function _init_particles(params::Parameters) 
    return rand(rng, 1:params.L, params.N)
end

function _init_lattice(particles, params::Parameters)
    lattice = zeros(eltype(particles), params.L)
    for x in particles
        lattice[x] += one(eltype(lattice))
    end
    return lattice
end

struct Chunks{T<:Integer}
    size::T
    num::T
    time_steps_between_saves::T
    function Chunks{T}(
        size::T, time_steps_between_saves::T, total_timesteps::T
    ) where {T<:Integer}
        num = div(total_timesteps, size*time_steps_between_saves)
        new(size, num, time_steps_between_saves)
    end
end
function Chunks(
    size::T, time_steps_between_saves::T, total_timesteps::T
) where T<:Integer
    return Chunks{T}(size, time_steps_between_saves, total_timesteps)
end