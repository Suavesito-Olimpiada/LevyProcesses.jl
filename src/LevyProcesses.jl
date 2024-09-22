module LevyProcesses

using Distributions
using LinearAlgebra
using InverseLaplace
using Integrals
using QuadGK
using NonlinearSolve
using StaticArrays

export sim

# indicator function of `cond`
𝟙(cond) = cond ? true : false
const ∑ = sum
const ∞ = Inf
arg(f, θ₀) = solve(NonlinearProblem((u, p) -> f(u), θ₀)).u
∫(f, a, b) = solve(IntegralProblem((u, p) -> f(u), (a, b)), QuadGKJL())
# inverse laplace transform, returns a function to evaluate
L⁻¹(f) = Weeks(f)


struct Emp{V<:AbstractVector}
    n::Int
    Z::V
    Emp(v::V) where {V} = new{V}(length(v), sort!(v))
end

# distrución empírica
(e::Emp)(x) = 1.0 / e.n * ∑(zᵢ -> 𝟙(zᵢ ≤ x), e.Z)
# (e::Emp)(x) = 1 / e.n * searchsortedlast(e.Z, x)


# needs q > 0
function sim(; c, σ, q, ϕf::Function, πf::Function, m, n)
    if q < 0 || false
        throw(ArgumentError("q must be positive, for the Gometric"))
    end
    Π(x) = ∫(πf, 0, x)
    ρ = arg(θ -> c * θ - ϕf(θ) - q + σ^2 * θ^2, 1)
    αᵨ = c + σ^2 * ρ
    Π̂(θ) = ∫(x -> (ℯ^(-θ * x) - ℯ^(-ρ * x)) * Π(x), 0, ∞)
    Ĥq(θ) = 1 / θ * Π̂(θ) / (ρ - θ) * 1 / (αᵨ - q / ρ)
    êσ(θ) = (αᵨ * σ^(-2)) / (αᵨ * σ^(-2) + θ) * 𝟙(σ > 0) + 𝟙(σ == 0)
    eσHq = L⁻¹(θ -> (êσ(θ) * Ĥq(θ))[])
    𝔻() = begin
        x = arg(x -> eσHq(x) - rand(Uniform(0, 1)), 1)
        x * 𝟙(0 < x < m)
    end
    μ = q / (ρ * αᵨ)
    Xq() = begin
        n = rand(Geometric(μ))
        n == 0 ? 0 : sum(𝔻() for _ in 1:n)
    end
    Yσρ() = rand(Exponential(αᵨ * σ^(-2) + ρ))
    Z() = Yσρ() + Xq()
    sample = [Z() for _ in 1:n]
    return ρ, αᵨ, Emp(sample)
end

end # module LevyProcesses
