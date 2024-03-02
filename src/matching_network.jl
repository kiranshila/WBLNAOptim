using CircuitNetworks, StaticArrays, LinearAlgebra

# Design of the matching network

#imn(L, ws, fs, tline_abcd_fn) = @tullio grad = Dual (*) s[j] := a2s <| tline_abcd_fn(fs[j], ws[i], L / length(ws))

function imn(L, ws::AbstractVector{T}, fs, tline_abcd_fn) where {T}
    NF = length(fs)
    NW = length(ws)
    δ = L / NW
    # Intermediate ABCD matrices
    eye = SMatrix{2,2,Complex{T}}(I)
    acc = fill(eye, NF)
    for j in eachindex(fs)
        for i in eachindex(ws)
            acc[j] *= tline_abcd_fn(fs[j], ws[i], δ)
        end
        acc[j] = a2s(acc[j])
    end
    acc
end

