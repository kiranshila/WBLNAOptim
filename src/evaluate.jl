# Functions to evaluate a result

const T₀ = 290.0f0;
gain2noise(gain) = T₀ * (1 - abs2(gain)) / abs2(gain)
nf2temp(nf) = T₀ * (nf - 1)
mag_s11_squared(s1, s2) = abs2(first(a2s(s2a(s1) * s2a(s2))))
mag_noise(s_match, nfmin, rn, sopt) = noise_figure(nfmin, rn, sopt, s_match[2, 2])
mag_noise(s_match, target) = mag_noise.(s_match, target.nfmin, target.rn, target.sopt)

mag_noise_Γout(Γout, nfmin, rn, sopt) = noise_figure(nfmin, rn, sopt, Γout)
mag_noise_Γout(Γout, target) = mag_noise_Γout.(Γout, target.nfmin, target.rn, target.sopt)