module HoppingRates

export hop_rates

hop_rates = Dict()

function simple_weak_strong_nr_hop_rate(species, n_A, n_B, alpha, chi)
    if species == 1
        return simple_weak_strong_nr_hop_rate_A(n_A, n_B, alpha, chi)
    elseif species == 2
        return simple_weak_strong_nr_hop_rate_B(n_A, n_B, alpha, chi)
    end
end

function simple_weak_strong_nr_hop_rate_A(n_A, n_B, alpha, chi)
    return n_B^(alpha) / ((n_A + n_B)^(chi))
end

function simple_weak_strong_nr_hop_rate_B(n_A, n_B, alpha, chi)
    return n_A^(-alpha) / ((n_A + n_B)^(chi))
end

function simple_weak_strong_nr_hop_rate_dt(N_A, N_B, alpha, chi)
    max_val = 0
    for nA in 1:N_A
        for nB in 1:N_B
	    hrA = nA*simple_weak_strong_nr_hop_rate_A(nA, nB, alpha, chi)
	    if hrA > max_val
                max_val = hrA
            end
	    hrB = nB*simple_weak_strong_nr_hop_rate_B(nA, nB, alpha, chi)
	    if hrB > max_val
                max_val = hrB
            end
        end
    end
    return 1/max_val
end

hop_rates["simple_weak_strong_nr_hop_rate"] = (
    simple_weak_strong_nr_hop_rate, simple_weak_strong_nr_hop_rate_dt
)

end
