"Construction of the Kravchuk matrix using three-term recurrence relation"

export kmatrix

begin
    k₁(x, k, j) = ((2*x)/√((2 + 2*j - k)*(k - 1)))
    k₂(k, j) = √(((3 + 2*j - k)*(k - 2))/((2 + 2*j - k)*(k - 1)))
end

function kmatrix(j)
    vk = zeros(Float64, (Integer(2*j + 1), Integer(2*j + 1)))
    Threads.@threads for x in -j:j
        vk[Integer(x + j + 1), 1] = 2.0^(-j) * √binomial(Integer(2*j), Integer(j + x))
        vk[Integer(x + j + 1), 2] = ((2*x)/√(2*j)) * vk[Integer(x + j + 1), 1]
        for k in 2:Integer(2*j+1)
            if k == 2
                vk[Integer(x + j + 1), k] = k₁(x, k, j) * vk[Integer(x + j + 1), k - 1]
            else
                vk[Integer(x + j + 1), k] = k₁(x, k, j) * vk[Integer(x + j + 1), k - 1] - k₂(k, j) * vk[Integer(x + j + 1), k - 2]
            end
        end
    end
    return vk
end