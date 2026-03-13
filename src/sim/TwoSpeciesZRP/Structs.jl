## ================================================== ##
##          Taking care of model parameters           ##
## ================================================== ##

abstract type Parameters end

struct OneSpeciesParameters{T<:Integer, S<:Real, U<:AbstractString} <: Parameters
    N::T
    L::T
    t_tot::T
    dt::S
    bc::U
end

struct TwoSpeciesParameters{T<:Integer, S<:Real, U<:AbstractString} <: Parameters
    N_A::T
    N_B::T
    N::T
    L::T
    t_tot::T
    dt::S
    alpha::S
    chi::S
    bc::U
end

abstract type State end

struct TwoSpeciesState{T<:Integer} <: State
    particles_A::Vector{T}
    particles_B::Vector{T}
    lattice_A::Vector{T}
    lattice_B::Vector{T}
end

function TwoSpeciesState(prm::TwoSpeciesParameters)
    particles_A = _init_particles(prm.N_A, prm.L)
    particles_B = _init_particles(prm.N_B, prm.L)

    lattice_A = _init_lattice(particles_A, prm.L)
    lattice_B = _init_lattice(particles_B, prm.L)

    return TwoSpeciesState(particles_A, particles_B, lattice_A, lattice_B)
end

struct ModelFunctions{F, G <: Function}
    hop_rate::F
    bc::G
end

function _init_particles(num_particles, lattice_size) 
    return rand(rng, 1:lattice_size, num_particles)
end

function _init_lattice(particles, lattice_size)
    lattice = zeros(eltype(particles), lattice_size)
    for x in particles
        lattice[x] += one(eltype(lattice))
    end
    return lattice
end

# TODO: Maybe move this to filehandling?
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