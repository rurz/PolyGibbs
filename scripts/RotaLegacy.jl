using Pkg
Pkg.activate(Base.current_project())

cd("..")

# This environment is for development #############################################################
#begin
#    using PyPlot
#    using LaTeXStrings
#    #matplotlib.use("TkAgg") #To use in Pluto.jl
#    rc("text", usetex = true)
#    ion()
#	pygui(true)
# end
###################################################################################################

using DelimitedFiles

const β = π/2
const θ = π/8

@info "Retrieving the two-dimensional field"
input = readdlm("data/origin/letterR/R_50.dat")

j = (size(input)[1] - 1)/2;

N = Integer(2*j + 1);

@info "The size of the field is $N"

#include("Krav.jl");

#K = reshape([Φ(n, q, j) for n in 0:2 * j for q in -j:j], Integer(2*j + 1), Integer(2*j + 1));

begin
    k₁(x, k, j) = ((2*x)/√((2 + 2*j - k)*(k - 1)))
    k₂(x, k, j) = √(((3 + 2*j - k)*(k - 2))/((2 + 2*j - k)*(k - 1)))
end;

vk = zeros(Float64, (Integer(2*j + 1), Integer(2*j + 1)));

Threads.@threads for x in -j:j
    vk[Integer(x + j + 1), 1] = 2^(-j) * √binomial(Integer(2*j), Integer(j + x))
    vk[Integer(x + j + 1), 2] = ((2*x)/√(2*j)) * vk[Integer(x + j + 1), 1]
    for k in 2:Integer(2*j+1)
        if k == 2
            vk[Integer(x + j + 1), k] = k₁(x, k, j) * vk[Integer(x + j + 1), k - 1]
        else
            vk[Integer(x + j + 1), k] = k₁(x, k, j) * vk[Integer(x + j + 1), k - 1] - k₂(x, k, j) * vk[Integer(x + j + 1), k - 2]
        end
    end
end

K = vk;

@info "Performing rotation at θ = $θ"

dmodes = zeros(ComplexF64, Integer(2 * j + 1), Integer(2 * j + 1));

begin
    coeff₁(n, k) = (2 * (-n/2 - (k - n/2 - 1) * cos(β)) * csc(β)) / √(k * (n - k + 1))

    coeff₂(n, k) = √((k - 1) * (n + 2 - k) / (k * (n + 1 - k)))

    coeff₃(n, k, l) = (2 * (1 - k + n/2 + (l - n/2 - 1) * cos(β)) * csc(β)) / √(l * (n - l + 1))

    coeff₄(n, k) = √((k - 1) * (n + 2 - k) / (k * (n + 1 - k)))
end;

Threads.@threads for n in 0:Integer(2*j)
    for m in -n:2:n
        vdd = zeros(Float64, (Integer(2 * j + 1), Integer(2 * j + 1)))
        vdd[1, 1] = cos(β/2)^n
        vdd[2, 1] = √n * (cos(β) - 1) * csc(β) * vdd[1, 1]
        for k in 2:n
            vdd[k + 1, 1] = coeff₁(n, k) * vdd[k, 1] - coeff₂(n ,k) * vdd[k - 1, 1]
        end
        for k in 1:(n + 1)
            for l in 1:n
                if l == 1
                    vdd[k, l + 1] = coeff₃(n, k, l) * vdd[k, l]
                else
                    vdd[k, l + 1] = coeff₃(n, k, l) * vdd[k, l] - coeff₄(n, l) * vdd[k, l - 1]
                end
            end
        end
        vld = zeros(ComplexF64, (Integer(2 * j + 1), Integer(2 * j + 1)))
        for x in 1:Integer(2*j+1)
            for y in 1:Integer(2*j+1)
                vld[x, y] = (-1)^((abs(m) - m)/2) * sum([K[x, nx + 1] * K[y, n - nx + 1] * (-1im)^(n - nx) * vdd[nx + 1, Integer((m + n)/2 + 1)] for nx in 0:n])
            end
        end
        vlD = zeros(ComplexF64, (Integer(2 * j + 1), Integer(2 * j + 1)))
        for x in 1:Integer(2*j+1)
            for y in 1:Integer(2*j+1)
                vlD[x, y] = conj(vld[x, y]) * sum([exp(-1im * m * θ) * vld[q, p] * input[q, p] for p in 1:Integer(2*j+1) for q in 1:Integer(2*j+1)])
            end
        end
        for x in 1:Integer(2*j+1)
            for y in 1:Integer(2*j+1)
                dmodes[x, y] = dmodes[x, y] + vlD[x, y]
            end
        end
    end
end

umodes = zeros(ComplexF64, Integer(2 * j + 1), Integer(2 * j + 1));

begin
    coeff₅(n, k, j) = (2 * ((1/2) * (n - 4 * j) - (k - 1 + (1/2) * (n - 4 * j)) * cos(β)) * csc(β)) / √(k * (4 * j + 1 - k - n))

    coeff₆(n, k, j) = √(((k - 1) * (4 * j + 2 - k - n)) / (k * (4 * j + 1 - k - n)))

    coeff₇(n, k, l, j) = (2 * (1 - k + (1/2) * (4 * j - n) + (l - 1 + (1/2) * (n - 4 * j)) * cos(β)) * csc(β)) / √(l * (4 * j + 1 - l - n))

    coeff₈(n, k, j) = √(((k - 1) * (4 * j + 2 - k - n)) / (k * (4 * j + 1 - k - n)))
end

Threads.@threads for n in Integer(4*j):-1:Integer(2*j)
    for m in Integer(-(4*j - n)):2:Integer(4 * j - n)
        vdu = zeros(Float64, (Integer(2 * j + 1), Integer(2 * j + 1)))
        vdu[1, 1] = cos(β/2)^(4 * j - n)
        vdu[2, 1] = √(4 * j - n) * (cos(β) - 1) * csc(β) * vdu[1, 1]
        for k in 2:Integer((4 * j - n))
            vdu[Integer(k + 1), Integer(1)] = coeff₅(n, k, j) * vdu[k, 1] - coeff₆(n, k, j) * vdu[k - 1, 1]
        end
        for k in 1:Integer((4 * j + 1 - n))
            for l in 1:Integer((4 * j - n))
                if l == 1
                    vdu[k, l + 1] = coeff₇(n, k, l, j) * vdu[k, l]
                else
                    vdu[k, l + 1] = coeff₇(n, k, l, j) * vdu[k, l] - coeff₈(n, l, j) * vdu[k, l - 1]
                end
            end
        end
        vlu = zeros(ComplexF64, (Integer(2 * j + 1), Integer(2 * j + 1)))
        for x in 1:Integer(2*j+1)
            for y in 1:Integer(2*j+1)
                vlu[x, y] = (-1)^((abs(m) - m)/2) * sum([K[x, Integer(2*j - nx + 1)] * K[y, Integer(nx - 2*j + n + 1)] * (-1im)^(n - nx) * vdu[nx + 1, Integer((m + 4*j - n)/2 + 1)] for nx in 0:Integer(4*j-n)])
            end
        end
        vlU = zeros(ComplexF64, (Integer(2 * j + 1), Integer(2 * j + 1)))
        for x in 1:Integer(2*j+1)
            for y in 1:Integer(2*j+1)
                vlU[x, y] = conj(vlu[x, y]) * sum([exp(-1im * m * θ) * vlu[q, p] * input[q, p] for p in 1:Integer(2*j+1) for q in 1:Integer(2*j+1)])
            end
        end
        for x in 1:Integer(2*j+1)
            for y in 1:Integer(2*j+1)
                umodes[x, y] = umodes[x, y] + vlU[x, y]
            end
        end
    end
end

imgrot = real(dmodes + umodes);

@info "Writting the results"
writedlm(joinpath(abspath(""), "data/target/letterR/R_50pi8.dat"), imgrot)
