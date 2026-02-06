#=
All basic structs are defined in this file including their 
setup functions.
=#

# Definition of parameters.
struct Parameters{T1<:Integer, T2<:Real, T3<:AbstractString}
    N::T1
    L::T1
    t::T1
    num_species::T1
    dt::T2
    bc::T3
end

# ZRP: Contains information about the particles positions, their species and
# and the occupancy of each lattice site. 
struct ZRP{T<:Integer}
    particles::Vector{T}
    species::Vector{T}
    lattice::Vector{T}
    spec_at_loc::Vector{T}
    current_t::T
end

function ZRP(prm::Parameters)
    particles = _init_particles(prm)    
    species = _init_species(prm)
    lattice = _init_lattice(particles, species, prm)
    specs_at_loc = zeros(eltype(lattice), prm.num_species)
    current_t = zero(eltype(particles))
    return ZRP(particles, species, lattice, specs_at_loc, current_t)
end

function _init_particles(prm)
    return rand(1:prm.L, prm.N)
end

function _init_species(prm)
    return rand(0:prm.num_species-1, prm.N)
end

function _init_lattice(particles, species, prm)
    lattice = zeros(eltype(particles), prm.L * prm.num_species)
    for (i, s) in zip(particles, species)
        lattice[i + s*prm.L] += one(eltype(lattice))
    end
    return lattice
end

# Chunk: Struct handling the chunking of data since it is more efficient to 
# write to HDF5 only every `chunk_t` steps. Additionally contains number of chunk
# `num_chunk`.
struct Chunk{T<:Integer} 
    particles::Matrix{T}
    lattice::Matrix{T}
    t::T
    num::T
end

function Chunk(zrp::ZRP, prm::Parameters, chunk_t::Integer) 
    chunk_particles = zeros(eltype(zrp.particles), (prm.N, chunk_t))
    chunk_particles[:, 1] = zrp.particles

    chunk_lattice = zeros(eltype(zrp.lattice), (prm.L * prm.num_species, chunk_t))
    chunk_lattice[:, 1] = zrp.lattice

    chunk_num = div(prm.t, chunk_t)

    return Chunk(chunk_particles, chunk_lattice, chunk_t, chunk_num)
    
end

# Improve computation by providing an extra structure:
struct Species{T<:Integer}
    nums::Vector{T}
end

function Species(prm::Parameters)
    nums = zeros(Int, prm.num_species)
    return Species(nums)
end


# Struct handling the boundary conditions
struct BC <: Function
    bc
end