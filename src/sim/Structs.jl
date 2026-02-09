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
    saving_time_step::T
end
