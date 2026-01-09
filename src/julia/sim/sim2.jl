
struct Parameters{T1<:Integer, T2<:Real, T3<:AbstractString}
    N::T1
    L::T1
    t::T1
    num_species::T1
    dt::T2
    bc::T3
end

struct ZRP{T<:Integer}
    particles::Vector{T}
    species::Vector{T}
    lattice::Vector{T}
end

function ZRP(prm::Parameters)
    particles = _init_particles(prm)    
    species = _init_species(prm)
    lattice = _init_lattice(particles, prm)
    return ZRP(particles, species, lattice)
end

function _init_particles(prm)
    return rand(1:prm.L, prm.N)
end

function _init_species(prm)
    return rand(0:prm.num_species-1, prm.N)
end

function _init_lattice(particles, prm)
    lattice = zeros(eltype(particles), prm.L)
    for i in particles
        lattice[i] += 1
    end
    return lattice
end