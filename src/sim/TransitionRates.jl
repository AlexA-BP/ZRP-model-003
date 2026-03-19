
module TransitionRates

using Random
include("TwoSpeciesZRP/Structs.jl")

export TransitionRatesDict

TransitionRatesDict = Dict()

function simple_weak_strong_nr_hop_rate(
    state, 
    model, 
    x, 
    spec,
)
    nA, nB = state.lattice[x, :]
    return simple_weak_strong_nr_hop_rate(spec, nA, nB, model.alpha, model.chi)
end


function simple_weak_strong_nr_hop_rate(species, n_A, n_B, alpha, chi)
    if species == 1
        return simple_weak_strong_nr_hop_rate_A(n_A, n_B, alpha, chi)
    elseif species == 2
        return simple_weak_strong_nr_hop_rate_B(n_A, n_B, alpha, chi)
    end
end

function simple_weak_strong_nr_hop_rate_A(n_A, n_B, alpha, chi)
    if n_A == 0
        return 0.
    end
    if n_B == 0
        return n_A^-(chi-1)
    end
    return n_A * n_B^(alpha) / ((n_A + n_B)^(chi))
end

function simple_weak_strong_nr_hop_rate_B(n_A, n_B, alpha, chi)
    if n_B == 0
        return 0.
    end
    if n_A == 0
        return n_B^-(chi-1)
    end
    return n_B * n_A^(-alpha) / ((n_A + n_B)^(chi))
end

function max_simple_weak_strong_nr_hop_rate(N_A, N_B, alpha, chi)
    max_val = 0
    for nA in 0:N_A, nB in 0:N_B
        if nA != 0
            hrA = simple_weak_strong_nr_hop_rate_A(nA, nB, alpha, chi)
            if hrA > max_val
                    max_val = hrA
            end
        end

        if nB != 0
            hrB = simple_weak_strong_nr_hop_rate_B(nA, nB, alpha, chi)
            if hrB > max_val
                    max_val = hrB
            end
        end
    end
    return max_val
end

function min_dt_simple_weak_strong_nr_hop_rate(N_A, N_B, alpha, chi)
    return 1/max_simple_weak_strong_nr_hop_rate(N_A, N_B, alpha, chi)
end

TransitionRatesDict["simple_weak_strong_nr_hop_rate"] = (
    simple_weak_strong_nr_hop_rate 
)

end
