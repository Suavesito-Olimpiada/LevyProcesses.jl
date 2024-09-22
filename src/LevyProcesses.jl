module LevyProcesses

using Distributions
using LinearAlgebra
using InverseLaplace
using Integrals
using QuadGK
using NonlinearSolve
using StaticArrays

export sim

𝟙(cond) = cond ? true : false
const ∑ = sum
const ∞ = Inf
arg(f, θ₀) = solve(NonlinearProblem((u, p) -> f(u), θ₀)).u
∫(f, a, b) = solve(IntegralProblem((u, p) -> f(u), (a, b)), QuadGKJL())
L⁻¹(f) = Weeks(f)


struct Emp{V<:AbstractVector}
    n::Int
    Z::V
    Emp(v::V) where {V} = new{V}(length(v), sort!(v))
end

(e::Emp)(x) = 1.0 / e.n * ∑(zᵢ -> 𝟙(zᵢ ≤ x), e.Z)
# (e::Emp)(x) = 1 / e.n * searchsortedlast(e.Z, x)


function sim(; c, σ, q, ϕf::Function, πf::Function, m, n)
    Π(x) = ∫(πf, 0, x)
    ρ = arg(ρ -> c * ρ - ϕf(ρ) - q + σ^2 * ρ^2, 1)
    @show ρ, c, σ
    αᵨ = c + σ^2 * ρ
    Π̂(θ, ρ) = ∫(x -> (ℯ^(-θ * x) - ℯ^(ρ * x)) * Π(x), 0, ∞)
    Ĥq(θ) = 1 / θ * Π̂(θ, ρ) / (ρ - θ) * 1 / (αᵨ - q / ρ)
    êσ(θ) = (αᵨ * σ^(-2)) / (αᵨ * σ^(-2) + θ) * 𝟙(σ > 0) + 𝟙(θ == 0)
    eσHq = L⁻¹(θ -> êσ(θ) * Ĥq(θ))
    𝔻() = begin
        x = arg(x -> eσHq(x) - rand(Uniform(0, 1)), 1)
        x * 𝟙(0 < x < m)
    end
    μ = 1 / (ρ * αᵨ)
    Xq() = begin
        n = rand(Geometric(μ))
        n == 0 ? 0 : sum(𝔻() for _ in 1:n)
    end
    Yσρ() = rand(Exponential(αᵨ * σ^(-2) + ρ))
    Z() = Yσρ() + Xq()
    sample = [Z() for _ in 1:n]
    return Emp(sample)
end

end # module LevyProcesses
