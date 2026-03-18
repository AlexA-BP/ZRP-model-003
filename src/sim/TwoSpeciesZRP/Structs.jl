## ================================================== ##
##          Taking care of model parameters           ##
## ================================================== ##

abstract type AbstractParameters end

struct ModelParameters{T<:Integer, S<:Real, U<:AbstractString} <: AbstractParameters
    Ns::Vector{T}
    num_spec::T
    N::T
    L::T
    phys_t::S
    num_tsteps::T
    dt::S
    bc::U
    alpha::S
    chi::S
    hop_directions::Tuple{T, T}
end


struct ModelFunctions{F1, F2 <: Function} <: AbstractParameters
    hop_rate::F1
    bc::F2
end

abstract type AbstractState end

struct NSpeciesState{T <: Integer} <: AbstractState
    particles::Vector{T}
    species::Vector{T}
    lattice::Matrix{T}
end

function NSpeciesState(model::ModelParameters, rng::AbstractRNG)
    particles = _init_particles(model.N, model.L, rng)
    species = _init_species(particles, model.Ns)
    lattice = _init_lattice(particles, species, model.L)
    return NSpeciesState(particles, species, lattice)
end

function _init_particles(num_particles, lattice_size, rng::AbstractRNG) 
    return rand(rng, 1:lattice_size, num_particles)
end

function _init_species(particles, nums_species)
    species = Vector{Int}(undef, size(particles))
    spec_names = collect(1:length(nums_species))
    for k in eachindex(species)
        i = pidgeon_hole_element_into_array(k, cumsum(nums_species))    
        species[k] = spec_names[i] 
    end
    return species
end

"""
For a sorted array of size n: a_1 <= a_2 <= ... <= a_n, pidgeon hole an element 
x into it. Returns the index i, s.t. a_i <= x <= a_{i+1}. Returns the last index 
n, if x > a_n. 
"""
function pidgeon_hole_element_into_array(x, arr)::Int
    for i in eachindex(arr)
        if x <= arr[i]
            return i
        end
    end
    return length(arr)
end

function _init_lattice(particles, species, L)
    lattice = zeros(eltype(particles), (L, 2))
    for (x, s) in zip(particles, species)
        lattice[x, s] += 1
    end        
    return lattice
end

# TODO: Maybe move this to filehandling?
struct Chunks{T<:Integer}
    size::T
    num::T
    time_steps_between_saves::T
    function Chunks{T}(
        size::T, time_steps_between_saves::T, total_timesteps::T, dt::S
    ) where {T<:Integer, S<:Real}
        num = div(total_timesteps, size*time_steps_between_saves)
        time_steps_between_saves = floor(time_steps_between_saves / dt)
        new(size, num, time_steps_between_saves)
    end
end
function Chunks(
    size::T, time_steps_between_saves::T, total_timesteps::T, dt::S
) where {T<:Integer, S<:Real}
    return Chunks{T}(size, time_steps_between_saves, total_timesteps, dt)
end

# 