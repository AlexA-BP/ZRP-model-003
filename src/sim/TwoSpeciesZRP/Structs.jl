

abstract type AbstractParameters end

struct ModelParameters{T<:Integer, S<:Real, U<:AbstractString} <: AbstractParameters
    # supplied directly
    Ns::Vector{T}
    L::T
    phys_t::S
    # determined by model
    num_spec::T
    N::T
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

mutable struct NSpeciesState{T <: Integer, S <: Real} <: AbstractState
    particles::Vector{T}
    species::Vector{T}
    lattice::Matrix{T}
    t::S
    dt::S
    trans_rates::Matrix{S}
    r0::S
end

function NSpeciesState(
    model::ModelParameters, 
    modelfunc::ModelFunctions, 
    rng::AbstractRNG
)
    particles = _init_particles(model.N, model.L, rng)
    species = _init_species(particles, model.Ns)
    lattice = _init_lattice(particles, species, model.L)
    trans_rates = _init_trans_rates(
        particles, species, lattice, model, modelfunc
    )
    r0 = maximum(trans_rates)
    return NSpeciesState(particles, species, lattice, 0., 0., trans_rates, r0)
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

function _init_trans_rates(
    particles, 
    species, 
    lattice, 
    model::ModelParameters, 
    modelfunc::ModelFunctions
)
    trans_rates = Matrix{Float64}(undef, size(lattice))
    for i in 1:model.L
        nA = lattice[i, 1] 
        nB = lattice[i, 2]
        
        trans_rates[i, 1] = modelfunc.hop_rate(
            1, nA, nB, model.alpha, model.chi
        )

        trans_rates[i, 2] = modelfunc.hop_rate(
            2, nA, nB, model.alpha, model.chi
        )
    end
    return trans_rates
end

# TODO: Maybe move this to filehandling?

abstract type AbstractChunks end

struct NaiveChunks{T<:Integer} <: AbstractChunks
    size::T
    num::T
    snapshot_phys_time_diff::T
    function NaiveChunks{T}(
        size::T, snapshot_phys_time_diff::T, total_timesteps::T, dt::S
    ) where {T<:Integer, S<:Real}
        num = div(total_timesteps, size*snapshot_phys_time_diff)
        snapshot_phys_time_diff = floor(snapshot_phys_time_diff / dt)
        new(size, num, snapshot_phys_time_diff)
    end
end
function NaiveChunks(
    size::T, snapshot_phys_time_diff::T, total_timesteps::T, dt::S
) where {T<:Integer, S<:Real}
    return Chunks{T}(size, snapshot_phys_time_diff, total_timesteps, dt)
end

struct Chunks{T<:Integer} <: AbstractChunks
    size::T
    num::T
    snapshot_phys_time_diff::T
end