using Revise

function test(N, L, t)
    v = rand(1:L, N)
    A = zeros(eltype(v), (N, t))
    
    update1!(A, v, N)

    return nothing
end


function update1!(A, v, N, t)
    for j in 1:t
        A[:, j] = v
    end
    return nothing
end

function update2!(A, v, N, t)
    for j in 1:t
        for i in 1:N
            A[i, j] = v[i]
        end
    end
    return nothing
end
