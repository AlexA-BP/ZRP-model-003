function update!(zrp::ZRP, prm::Parameters)
    for _ in 1:prm.N
        x = zrp.particles[rand(1:prm.N)]
        n = zrp.lattice[x]

        p = hop_rate(n)
    end
end

function hop_rate(n) end