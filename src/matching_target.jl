# Matching Target tline_model

using Tullio

struct TwoPortNetwork{T,S<:AbstractMatrix{Complex{T}}}
    freqs::Vector{T}
    s::Vector{S}
    nfmin::Vector{T}
    rn::Vector{T}
    sopt::Vector{Complex{T}}
end

Base.show(io::IO, _::TwoPortNetwork) = show(io, "A 2-Port Network")

function parse_ads_complex(s)
    mag, ang = split(s, "/")
    mag = strip(mag)
    ang = strip(ang)
    mag = parse(Float32, mag)
    ang = parse(Float32, ang)
    mag * cispi(ang / 180.0f0)
end

function TwoPortNetwork(csv_filename)
    file = CSV.File(csv_filename; comment="#")
    freqs = Float32.(file["freq"])
    nfmin = @. Float32(10^(file["nfmin"] / 10))
    rn = Float32.(file["rn"])
    sopt = file["sopt"] .|> parse_ads_complex
    @tullio S[i] := @SMatrix [
        parse_ads_complex(file[i].s11) parse_ads_complex(file[i].s12);
        parse_ads_complex(file[i].s21) parse_ads_complex(file[i].s22)
    ]
    TwoPortNetwork(freqs, S, nfmin, rn, sopt)
end

const TARGET = TwoPortNetwork("$(@__DIR__)/../data/LNA_Noise_S.csv")
const DIE_4F50 = TwoPortNetwork("$(@__DIR__)/../data/4F50_Noise_S.csv")