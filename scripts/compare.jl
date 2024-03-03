using CSV

hfss = CSV.File("$(@__DIR__)/../simulations/HFSS_optim10.csv"; types=[Float64, ComplexF64, ComplexF64, ComplexF64, ComplexF64])
fs = range(0.7e9, 2e9, 201)

hfss_s = [@SMatrix [r[1].s11 r[1].s12; r[1].s21 r[1].s22] for r in eachrow(hfss)]

noise_temp = nf2temp.(mag_noise(hfss_s, TARGET))
insertion_noise = @. gain2noise(available_gain(hfss_s, 0.0))
total_noise = @. insertion_noise + noise_temp
input_tmin = @. insertion_noise + nf2temp(TARGET.nfmin)
input_rl = 10 .* log10.(mag_s11_squared.(hfss_s, TARGET.s))
plot_freqs = TARGET.freqs .* 1e-9
plot(plot_freqs, total_noise, label="Noise")
plot!(plot_freqs, input_tmin, label="Minimum Noise", ylims=(0, 20), xlabel="Freq (GHz)", ylabel="Noise (K)", legend=false)

θ = readdlm("$(@__DIR__)/../shape_results/optim_10.csv", ',')
s_match = imn(θ[1], θ[2:end], TARGET.freqs, em_abcd)
noise_temp = nf2temp.(mag_noise(s_match, TARGET))
insertion_noise = @. gain2noise(available_gain(s_match, 0.0))
total_noise = @. insertion_noise + noise_temp
input_tmin = @. insertion_noise + nf2temp(TARGET.nfmin)
input_rl = 10 .* log10.(mag_s11_squared.(s_match, TARGET.s))
plot_freqs = TARGET.freqs .* 1e-9
plot!(plot_freqs, total_noise, label="Noise")
plot!(plot_freqs, input_tmin, label="Minimum Noise", ylims=(0, 20), xlabel="Freq (GHz)", ylabel="Noise (K)", legend=false)