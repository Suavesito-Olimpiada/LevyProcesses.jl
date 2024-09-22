using LinearAlgebra
using LevyProcesses
using ForwardDiff
using MittagLeffler
using NonlinearSolve
using LaTeXStrings
using Plots


const LP = LevyProcesses

# μ = 2.0
# σ = 1.0
# q = 0.0  # usually q > 0, at least q ≥ 0
# # this does not correspond to a "normal" case
# ϕf(θ) = LP.𝟙(θ > 0) * (θ^2*σ^2/2.0 + μ*θ)
# πf(x) = 0.0

# μ = 2.0
# σ = 1.0
# q = 0.0  # usually q > 0, at least q ≥ 0
# # this correspond to the original case by Ehyter
# ϕf(x) = LP.𝟙(x > 0) * exp(-x)
# πf(x) = 1.0 - 1.0 / (1.0 + x)
# wq(x) = 2/(√(2*q*σ^2+μ)) * ℯ^(-μ*x/σ^2) * sinh((x/σ^2)*√(2*q*σ^2+μ)) * (x > 0)

# #  this needs that `c-λ/μ > 0`
# c = 1.0  # drift
# λ = 1.0  # rate of arrival
# μ = 2.0  # exponentially distributed jumps
# # usually q > 0, at least q ≥ 0. this is where the q of q-scale, and, therefore, the q W(q)
# # comes from
# q = 0.0
# # the characteristic triplet (σ, ϕf, πf)
# σ = 0.0
# ϕf(θ) = LP.𝟙(θ > 0) * (c * θ - λ * (1 - μ / (μ + θ)))
# # ρ = arg(θ -> c * θ - ϕf(θ) - q + σ^2 * θ^2, 1)
# πf(x) = exp(-μ * x)
# wq(x) = 1 / c * (1 + λ / (c * μ - λ) * (1 - ℯ^(μ - c^-1 * λ) * x))

# # case from Ronnie Loeffen's Stochastic control for spectrally negative Levy processes)
# # ch 4.10, pp 73.
# N = 10
# Aₖ = normalize!(rand(N))
# αₖ = 1
# # μ = 2.0
# # σ = 1.0
# # q = 0.0  # usually q > 0, at least q ≥ 0
# # # this does not correspond to a "normal" case
# # ϕf(θ) = LP.𝟙(θ > 0) * (θ^2*σ^2/2.0 + μ*θ)
# # πf(x) = 0.0

# # case from Ronnie Loeffen's Stochastic control for spectrally negative Levy processes)
# # ch 1.6, pp 73.
# α = 1  # > 0
# λ = 10
# c = μ = 21.4
# σ₁ = 1.4
# σ₂ = 2
# q = 0.1  # usually q > 0, at least q ≥ 0
# # ψ(u) in the tesis
# _ϕf(θ) = c*θ + θ^2*σ^2/2.0 + λ*(α/(θ+α)^2 - 1)
# ϕf(θ) = LP.𝟙(θ > 0) * _ϕf(θ)
# πf(x) = λ*α^2*x*exp(-α*x)
#
# f(θ) = (_ϕf(θ) - q)*(α + θ)^2
# θⱼ = solve(NonlinearProblem((u, p) -> f(u), 1))


m = 100


for n in (10, 100, 1000)
    @time ρ, αᵨ, emp = sim(; c, σ, q, ϕf, πf, m, n)

    x = -0.1:0.1:10
    y = emp.(x)
    Y = @. (ρ * αᵨ * σ^(-2)) / (q * (αᵨ * σ^(-2) + ρ)) * ℯ^(ρ * x) * y
    Yw = wq.(x)

    plt = plot(x, y)
    plt2 = plot(x, Y)

    savefig(plt, "emp_known_p$n.png")
    savefig(plt2, "emp_known_$n.png")
end

#  1.969966 seconds (4.90 M allocations: 159.307 MiB, 2.11% gc time)
# 11.379837 seconds (24.47 M allocations: 492.415 MiB, 0.84% gc time)
# 89.011660 seconds (192.06 M allocations: 3.266 GiB, 0.78% gc time)
