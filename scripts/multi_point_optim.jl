using CircuitNetworks, Statistics, Enzyme, ForwardDiff, Optimization, OptimizationOptimJL

include("$(@__DIR__)/../src/tline_model.jl")
include("$(@__DIR__)/../src/matching_network.jl")
include("$(@__DIR__)/../src/matching_target.jl")
include("$(@__DIR__)/../src/evaluate.jl")
include("$(@__DIR__)/../src/plotting.jl")

function cost(L, widths, target)
    s_match = imn(L, widths, target.freqs, em_abcd)
    noise_temp = nf2temp.(mag_noise(s_match, target))
    insertion_noise = @. gain2noise(available_gain(s_match, 0.0))
    total_noise = insertion_noise .+ noise_temp
    input_rl = mag_s11_squared.(s_match, target.s)

    rl_cost = mean(input_rl)
    noise_cost = (mean(total_noise) / 8.5 - 1)^2

    100 * noise_cost + rl_cost
end

cost(θ, target) = cost(θ[1], θ[2:end], target)

callback = function (state, cost)
    θ = state.u
    @show cost
    display(plot_performance(θ[1], θ[2:end], TARGET))
    return false
end

# Solving
N = 100
θ₀ = [90e-3, fill(5e-3, N)...]
lb = [20e-3, fill(0.2e-3, N)...]
ub = [120e-3, fill(20e-3, N)...]

optf = OptimizationFunction(cost, Optimization.AutoEnzyme())
prob = OptimizationProblem(optf, θ₀, TARGET; callback=callback, lb=lb, ub=ub)
sol = solve(prob, BFGS())

cost(θ) = cost(θ, TARGET)