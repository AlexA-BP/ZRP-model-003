module HoppingRates

export hop_rates

hop_rates = Dict()

function simple_weak_strong_nr_hop_rate(species, n_A, n_B, alpha, chi)
    if species == "A"
        return simple_weak_strong_nr_hop_rate_A(n_A, n_B, alpha, chi)
    elseif species == "B"
        return simple_weak_strong_nr_hop_rate_B(n_A, n_B, alpha, chi)
    end
end

function simple_weak_strong_nr_hop_rate_A(n_A, n_B, alpha, chi)
    return n_B^(alpha) / ((n_A + n_B)^(chi))
end

function simple_weak_strong_nr_hop_rate_B(n_A, n_B, alpha, chi)
    return n_A^(-alpha) / ((n_A + n_B)^(chi))
end

hop_rates["simple_weak_strong_nr_hop_rate"] = simple_weak_strong_nr_hop_rate
