using CircuitNetworks, Statistics, Optimization, OptimizationOptimJL, CSV, DelimitedFiles
using OptimizationNLopt
using ForwardDiff, PolyesterForwardDiff

include("$(@__DIR__)/../src/tline_model.jl")
include("$(@__DIR__)/../src/coax_model.jl")
include("$(@__DIR__)/../src/matching_network.jl")
include("$(@__DIR__)/../src/matching_target.jl")
include("$(@__DIR__)/../src/evaluate.jl")
include("$(@__DIR__)/../src/plotting.jl")

function cost(L, inner_diameters, outer_diameter, target)
    s_match = coax_imn(L, inner_diameters ./ 2, outer_diameter / 2, target.freqs)
    input_stub_abcds = em_abcd.(target.freqs, 1.5e-3, 1.5e-3)
    s_match = a2s.(s2a.(s_match) .* input_stub_abcds)
    noise_temp = nf2temp.(mag_noise(s_match, target))
    insertion_noise = @. gain2noise(available_gain(s_match, 0.0))
    total_noise = noise_temp .+ insertion_noise
    total_minimum_noise = nf2temp.(target.nfmin) .+ insertion_noise
    input_rl = mag_s11_squared.(s_match, target.s)
    rl_cost = mean(input_rl)
    noise_cost = (mean(total_noise) / 5 - 1)^2
    noise_mismatch_cost = mean(@. (total_noise / total_minimum_noise - 1)^2)
    10 * noise_cost + 3 * noise_mismatch_cost + 1 * rl_cost

    # Regularize smoothness
    #step_cost = mean(abs2.(widths[2:end] .- widths[1:end-1]))

    #noise_mismatch_cost + 5000 * step_cost
end

cost(θ, target) = cost(θ[1], θ[2:end], 20e-3, target)

callback = function (state, cost)
    θ = state.u
    @show cost
    display(plot_performance_coax(θ[1], θ[2:end], 20e-3, TARGET))
    #cost < 0.135
    return false
end

# Solving
N = 100
θ₀ = [90e-3, fill(2e-3, N)...]
lb = [20e-3, fill(0.5e-3, N)...]
ub = [120e-3, fill(10e-3, N)...]
optf = OptimizationFunction(cost, Optimization.AutoForwardDiff())
prob = OptimizationProblem(optf, θ₀, TARGET; callback=callback, lb=lb, ub=ub)
sol = solve(prob, NLopt.LD_LBFGS())

writedlm("$(@__DIR__)/../shape_results/optim_$(N).csv", sol.u, ',')
savefig(plot_performance(sol.u[1], sol.u[2:end], TARGET), "$(@__DIR__)/../plots/optim_$(N).png")

