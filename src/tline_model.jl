# Create model interpolated from HFSS 2D EM simulation data
using Interpolations, CSV, CircuitNetworks, StaticArrays

function abcd_tline(γ, z₀, d)
    sh = sinh(γ * d)
    ch = cosh(γ * d)
    @SMatrix [ch sh*z₀; sh/z₀ ch]
end

function abcd_tline(r, l, g, c, f, d)
    gamma = CircuitNetworks.γ(r, l, g, c, f)
    z0 = CircuitNetworks.Z₀(r, l, g, c, f)
    abcd_tline(gamma, z0, d)
end

interp(ys, xs) = scale(interpolate(ys, BSpline(Linear())), xs)

rlgc_data = CSV.File("$(@__DIR__)/../data/2D_EM_RLGC.csv"; comment="#")
rlgc_f_range = range(0.6e9, 2.1e9, 51)
rlgc_w_range = range(0.2e-3, 15e-3, step=0.1e-3)
data_range = (rlgc_f_range, rlgc_w_range)
data_shape = (length(rlgc_f_range), length(rlgc_w_range))

# Extract and rescale data
rlgc_r = reshape(rlgc_data["r_ohm"], data_shape) .|> Float32
rlgc_l = reshape(rlgc_data["l_uH"], data_shape) .* 1e-6 .|> Float32
rlgc_g = reshape(rlgc_data["g_uS"], data_shape) .* 1e-6 .|> Float32
rlgc_c = reshape(rlgc_data["c_pF"], data_shape) .* 1e-12 .|> Float32

const r_model = interp(rlgc_r, data_range)
const l_model = interp(rlgc_l, data_range)
const g_model = interp(rlgc_g, data_range)
const c_model = interp(rlgc_c, data_range)

function em_abcd(f, w, d)
    r = r_model(f, w)
    l = l_model(f, w)
    g = g_model(f, w)
    c = c_model(f, w)
    abcd_tline(r, l, g, c, f, d)
end