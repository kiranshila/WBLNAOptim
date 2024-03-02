# Plotting routines for evaluation

using Plots

default(:dpi, 600)
default(:size, (800, 500))

function plot_performance(L, widths, target)
    s_match = imn(L, widths, target.freqs, em_abcd)
    noise_temp = nf2temp.(mag_noise(s_match, target))
    insertion_noise = @. gain2noise(available_gain(s_match, 0.0))
    total_noise = @. insertion_noise + noise_temp
    input_tmin = @. insertion_noise + nf2temp(target.nfmin)
    input_rl = 10 .* log10.(mag_s11_squared.(s_match, target.s))
    plot_freqs = target.freqs .* 1e-9
    plot(plot_freqs, total_noise, label="Noise")
    noise_plot = plot!(plot_freqs, input_tmin, label="Minimum Noise", ylims=(0, 20), xlabel="Freq (GHz)", ylabel="Noise (K)", legend=false)
    s11_plot = plot(plot_freqs, input_rl; label="S11", ylims=(-20, 0), xlabel="Freq (GHz)", ylabel="S11 (dB)", legend=false)
    points = range(0, L, length(widths)) .* 1000
    plot(points, widths ./ 2 .* 1000, aspect_ratio=1, linetype=:steppre)
    shape_plot = plot!(points, -widths ./ 2 .* 1000, aspect_ratio=1, xlims=(0, 120), ylims=(-20, 20), legend=false, linetype=:steppre)
    l = @layout [
        [a{0.5w} b{0.5w}]
        c{1.0w}
    ]
    plot(noise_plot, s11_plot, shape_plot, layout=l)
end