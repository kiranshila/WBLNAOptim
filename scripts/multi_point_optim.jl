using CircuitNetworks, Statistics, Optimization, OptimizationOptimJL, CSV, DelimitedFiles
using OptimizationNLopt
using ForwardDiff, PolyesterForwardDiff, Enzyme

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
    total_minimum_noise = insertion_noise .+ nf2temp.(target.nfmin)
    input_rl = mag_s11_squared.(s_match, target.s)
    rl_cost = mean(input_rl)
    noise_cost = (mean(total_noise) / 5 - 1)^2
    noise_mismatch_cost = mean(@. (total_noise / total_minimum_noise - 1)^2)
    5 * noise_cost + 3 * noise_mismatch_cost + 10 * rl_cost
end

cost(θ, target) = cost(θ[1], θ[2:end], target)

callback = function (state, cost)
    θ = state.u
    @show cost
    display(plot_performance(θ[1], θ[2:end], TARGET))
    #cost < 0.135
    return false
end

# Solving
for N in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 25, 50, 100]
    θ₀ = [90e-3, fill(3e-3, N)...]
    lb = [20e-3, fill(0.2e-3, N)...]
    ub = [120e-3, fill(20e-3, N)...]
    optf = OptimizationFunction(cost, Optimization.AutoForwardDiff())
    prob = OptimizationProblem(optf, θ₀, TARGET; callback=callback, lb=lb, ub=ub)
    sol = solve(prob, NLopt.LD_LBFGS())
    writedlm("$(@__DIR__)/../shape_results/optim_$(N).csv", sol.u, ',')
    savefig(plot_performance(sol.u[1], sol.u[2:end], TARGET), "$(@__DIR__)/../plots/optim_$(N).png")
end