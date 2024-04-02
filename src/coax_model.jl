const ε₀ = 8.85418782e-12
const μ₀ = 1.25663706e-6

σ_silver = 6.3e7
σ_copper = 5.8e7
σ_aluminum = 3.5e7

μᵣ_silver = 0.99998

function coax_rlgc(a, b, μ, ε, σ, μc, σc, freq)
    Rs = sqrt(π * freq * μc / σc)
    lba = log(b / a)
    R = Rs / (2π) * (1 / a + 1 / b)
    L = μ / (2π) * lba
    G = 2π * σ / lba
    C = (2π * ε) / lba
    (R, L, G, C)
end

silver_coax_rlgc(a, b, freq) = coax_rlgc(a, b, μ₀, ε₀, 0, μᵣ_silver * μ₀, σ_silver, freq)

silver_coax_abcd(f, a, b, l) = abcd_tline(silver_coax_rlgc(a, b, f)..., f, l)

function coax_imn(L, as::AbstractVector{T}, b, fs) where {T}
    NF = length(fs)
    NA = length(as)
    δ = L / NA
    # Intermediate ABCD matrices
    eye = SMatrix{2,2,Complex{T}}(I)
    acc = fill(eye, NF)
    for j in eachindex(fs)
        for i in eachindex(as)
            acc[j] *= silver_coax_abcd(fs[j], as[i], b, δ)
        end
        acc[j] = a2s(acc[j])
    end
    acc
end