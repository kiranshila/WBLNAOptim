using Enzyme, Zygote

f(x) = abs(sqrt(x + x * im))

Zygote.gradient(f, 1.0)