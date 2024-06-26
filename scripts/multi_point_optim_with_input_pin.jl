using CircuitNetworks, Statistics, Optimization, OptimizationOptimJL, CSV, DelimitedFiles
using OptimizationNLopt
using ForwardDiff, PolyesterForwardDiff

include("$(@__DIR__)/../src/tline_model.jl")
include("$(@__DIR__)/../src/matching_network.jl")
include("$(@__DIR__)/../src/matching_target.jl")
include("$(@__DIR__)/../src/evaluate.jl")
include("$(@__DIR__)/../src/plotting.jl")

function cost(L, widths, target)
    # Compute input matching network with N connector pin
    s_match = imn(L, widths, target.freqs, em_abcd)
    input_stub_abcds = em_abcd.(target.freqs, 1.5e-3, 1.5e-3)
    s_match = a2s.(s2a.(s_match) .* input_stub_abcds)

    # Compute noise from Γs and available gain of the input network
    noise_temp = noise_temperature.(mag_noise(s_match, target))
    insertion_noise = noise_temperature.(noise_figure.(s_match, 0.0, 290.0))

    # Add the noises to get the true noise and Tmin at the input
    total_noise = noise_temp .+ insertion_noise
    total_minimum_noise = noise_temperature.(target.nfmin) .+ insertion_noise

    # Compute return loss at the input
    input_rl = mag_s11_squared.(s_match, target.s)

    # Formulate the costs
    rl_cost = mean(input_rl)
    noise_cost = (mean(total_noise) / 5 - 1)^2
    noise_mismatch_cost = mean(@. (total_noise / total_minimum_noise - 1)^2)

    # Create the resulting cost to optimize
    10 * noise_cost + 3 * noise_mismatch_cost + 1 * rl_cost

    # Regularize smoothness
    #step_cost = mean(abs2.(widths[2:end] .- widths[1:end-1]))
    #10 * noise_cost + 3 * noise_mismatch_cost + 5000 * step_cost #+ 1 * rl_cost
end

cost(θ, target) = cost(θ[1], θ[2:end], target)

callback = function (state, cost)
    θ = state
    display(plot_performance(θ[1], θ[2:end], TARGET))
    #cost < 0.135
    return false
end

# Solving
N = 10
θ₀ = [90e-3, fill(2e-3, N)...]
lb = [20e-3, fill(0.2e-3, N)...]
ub = [200e-3, fill(15e-3, N)...]
optf = OptimizationFunction(cost, Optimization.AutoForwardDiff())
prob = OptimizationProblem(optf, θ₀, TARGET; callback=callback, lb=lb, ub=ub)
sol = solve(prob, NLopt.LD_LBFGS())

writedlm("$(@__DIR__)/../shape_results/optim_$(N).csv", sol.u, ',')
savefig(plot_performance(sol.u[1], sol.u[2:end], TARGET), "$(@__DIR__)/../plots/optim_$(N).png")

plot_performance(sol.u[1], sol.u[2:end], TARGET)