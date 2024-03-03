using Measurements, DelimitedFiles

θ = readdlm("$(@__DIR__)/../shape_results/optim_10.csv", ',')
L = θ[1]
ws = θ[2:end]

ws_uncertain = ws .± 0.0762e-3 # 3 mil width error
