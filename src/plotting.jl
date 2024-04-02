# Plotting routines for evaluation

using Plots

default(:dpi, 600)
default(:size, (800, 500))

function plot_performance(L, widths, target)
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

    plot_freqs = target.freqs .* 1e-9
    plot(plot_freqs, total_noise, label="Noise")
    plot!(plot_freqs, total_minimum_noise, label="Minimum Noise", ylims=(0, 20), xlabel="Freq (GHz)", ylabel="Noise (K)", legend=false)
    noise_plot = annotate!(1.2, 3, "Avg Noise: $(mean(total_noise))", 10)
    s11_plot = plot(plot_freqs, input_rl; label="S11", ylims=(-20, 0), xlabel="Freq (GHz)", ylabel="S11 (dB)", legend=false)

    points = range(0, L, length(widths) + 1) .* 1000
    include_end = [widths..., widths[end]]
    #plot(points, include_end .* 500, linetype=:steppost, color=:black)
    shape_plot = plot(points, -include_end .* 500,
        fillrange=include_end .* 500,
        aspect_ratio=1,
        xlims=(0, 120),
        ylims=(-20, 20),
        legend=false,
        linetype=:steppost,
        color=:grey)

    l = @layout [
        [a{0.5w} b{0.5w}]
        c{1.0w}
    ]
    plot(noise_plot, s11_plot, shape_plot, layout=l)
end

function plot_performance_100(L, Zfeed, widths, target)
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

    # FIXME
    input_rl = 10 .* log10.(mag_s11_squared.(s_match, target.s))

    plot_freqs = target.freqs .* 1e-9
    plot(plot_freqs, total_noise, label="Noise")
    plot!(plot_freqs, total_minimum_noise, label="Minimum Noise", ylims=(0, 20), xlabel="Freq (GHz)", ylabel="Noise (K)", legend=false)
    noise_plot = annotate!(1.2, 3, "Avg Noise: $(mean(total_noise))", 10)
    s11_plot = plot(plot_freqs, input_rl; label="S11", ylims=(-20, 0), xlabel="Freq (GHz)", ylabel="S11 (dB)", legend=false)

    points = range(0, L, length(widths) + 1) .* 1000
    include_end = [widths..., widths[end]]
    #plot(points, include_end .* 500, linetype=:steppost, color=:black)
    shape_plot = plot(points, -include_end .* 500,
        fillrange=include_end .* 500,
        aspect_ratio=1,
        xlims=(0, 120),
        ylims=(-20, 20),
        legend=false,
        linetype=:steppost,
        color=:grey)

    l = @layout [
        [a{0.5w} b{0.5w}]
        c{1.0w}
    ]
    plot(noise_plot, s11_plot, shape_plot, layout=l)
end

function plot_performance_coax(L, inner_diameters, outer_diameter, target)
    s_match = coax_imn(L, inner_diameters ./ 2, outer_diameter / 2, target.freqs)

    # Input pin
    input_stub_abcds = em_abcd.(target.freqs, 1.5e-3, 1.5e-3)
    s_match = a2s.(s2a.(s_match) .* input_stub_abcds)

    noise_temp = nf2temp.(mag_noise(s_match, target))
    insertion_noise = @. gain2noise(available_gain(s_match, 0.0))
    total_noise = @. insertion_noise + noise_temp
    input_tmin = @. insertion_noise + nf2temp(target.nfmin)
    input_rl = 10 .* log10.(mag_s11_squared.(s_match, target.s))
    plot_freqs = target.freqs .* 1e-9
    plot(plot_freqs, total_noise, label="Noise")
    plot!(plot_freqs, input_tmin, label="Minimum Noise", ylims=(0, 20), xlabel="Freq (GHz)", ylabel="Noise (K)", legend=false)
    noise_plot = annotate!(1.2, 3, "Avg Noise: $(mean(total_noise))", 10)
    s11_plot = plot(plot_freqs, input_rl; label="S11", ylims=(-20, 0), xlabel="Freq (GHz)", ylabel="S11 (dB)", legend=false)

    points = range(0, L, length(inner_diameters) + 1) .* 1000
    include_end = [inner_diameters..., inner_diameters[end]]
    #plot(points, include_end .* 500, linetype=:steppost, color=:black)
    shape_plot = plot(points, -include_end .* 500,
        fillrange=include_end .* 500,
        aspect_ratio=1,
        xlims=(0, 120),
        ylims=(-20, 20),
        legend=false,
        linetype=:steppost,
        color=:grey)

    l = @layout [
        [a{0.5w} b{0.5w}]
        c{1.0w}
    ]
    plot(noise_plot, s11_plot, shape_plot, layout=l)
end