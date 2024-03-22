using CircuitNetworks, Statistics, Optimization, OptimizationOptimJL, CSV, DelimitedFiles
using OptimizationNLopt
using ForwardDiff, PolyesterForwardDiff

include("$(@__DIR__)/../src/tline_model.jl")
include("$(@__DIR__)/../src/matching_network.jl")
include("$(@__DIR__)/../src/matching_target.jl")
include("$(@__DIR__)/../src/evaluate.jl")
include("$(@__DIR__)/../src/plotting.jl")

Γout(S, Γs) = S[2, 2] + (S[1, 2] * S[2, 1] * Γs) / (1 - S[1, 1] * Γs)

function cost(L, Zfeed, widths, target)
    # Noise
    s_match = imn(L, widths, target.freqs, em_abcd)
    input_stub_abcds = em_abcd.(target.freqs, 1.5e-3, 1.5e-3)
    s_match = a2s.(s2a.(s_match) .* input_stub_abcds)

    Γs = CircuitNetworks.Γ(Zfeed, 50)
    gamma_out = Γout.(s_match, Γs)
    noise_figure = mag_noise_Γout(gamma_out, target)
    noise_temp = nf2temp.(noise_figure)

    insertion_noise = @. gain2noise(available_gain(s_match, Γs))
    total_noise = insertion_noise .+ noise_temp
    total_minimum_noise = insertion_noise .+ nf2temp.(target.nfmin)

    noise_cost = (mean(total_noise) / 5 - 1)^2
    noise_mismatch_cost = mean(@. (total_noise / total_minimum_noise - 1)^2)

    # Return loss
    #input_rl = mag_s11_squared.(s_match, target.s)
    #rl_cost = mean(input_rl)

    # Total weighted cost
    10 * noise_cost + 3 * noise_mismatch_cost #+ 1 * rl_cost
end

cost(θ, target) = cost(θ[1], θ[2], θ[3:end], target)

callback = function (state, cost)
    θ = state.u
    @show cost
    display(plot_performance_100(θ[1], θ[2], θ[3:end], DIE_4F50))
    #cost < 0.135
    return false
end

# Solving
N = 25
θ₀ = [90e-3, 20, fill(1e-3, N)...]
lb = [20e-3, 19.0, fill(0.4e-3, N)...]
ub = [120e-3, 501.0, fill(15e-3, N)...]
optf = OptimizationFunction(cost, Optimization.AutoForwardDiff())
prob = OptimizationProblem(optf, θ₀, DIE_4F50; callback=callback, lb=lb, ub=ub)
sol = solve(prob, NLopt.LD_LBFGS())

#writedlm("$(@__DIR__)/../shape_results/optim_$(N).csv", sol.u, ',')
#savefig(plot_performance(sol.u[1], sol.u[2:end], TARGET), "$(@__DIR__)/../plots/optim_$(N).png")
