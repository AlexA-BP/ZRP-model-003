using StaticArrays

function update!(zrp::ZRP, prm::Parameters, bc::F) where F<:Function
    for _ in 1:prm.N
        k = rand(1:prm.N)
        xk = zrp.particles[k]
        sk = zrp.species[k]

        for s in 0:prm.num_species-1
            zrp.spec_at_loc[s+1] = zrp.lattice[xk + s*prm.L]
        end
        num_at_x = sum(zrp.spec_at_loc)

        hopprob = hop_rate(num_at_x)*prm.dt/num_at_x
        if rand() < hopprob
            hop!(zrp, prm, bc, k, xk, sk)            
        end
    end
    return nothing
end

function hop_rate(n)
    return 1. + 3.5/n
end

function hop_rate(ns::AbstractVector, n::Integer) 
    return 1.0/n
end

function hop!(zrp::ZRP, prm::Parameters, bc, i, x0, s)
    x1 = new_position(x0, prm, bc)
    hop_particles!(zrp, i, x1)    
    hop_lattice!(zrp, x0, x1, prm.L, s)
    return nothing
end

function new_position(x0, prm::Parameters, bc)
    x1 = x0 + rand((oneunit(x0), -oneunit(x0)))
    # x1 = x0 + oneunit(x0)
    x1 = bc(x1, prm.L)
    return x1    
end

function hop_particles!(zrp::ZRP, i, x1)
    zrp.particles[i] = x1
    return nothing
end

function hop_lattice!(zrp::ZRP, x0, x1, L, s)
    zrp.lattice[x0 + L*s] -= 1
    zrp.lattice[x1 + L*s] += 1
    return nothing
end

function apply_bc(x, prm::Parameters)
    if prm.bc in ["p", "periodic"]
        return apply_periodic_bc(x, prm.L)
    else
        error("Implementation error")
    end
end

function periodic_bc(x, L)
    if x == zero(x)
        return L
    elseif x == L+oneunit(x)
        return oneunit(x)
    end
    return x
end

function reflective_bc(x, L)
    if x == zero(x)
        return 2
    elseif x == L+oneunit(x)
        return L-oneunit(x)
    end
    return x
end