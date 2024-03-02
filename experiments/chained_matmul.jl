using CUDA, StaticArrays, BenchmarkTools

test_array_cpu = fill((@SMatrix rand(ComplexF32, 2, 2)), 1000)
test_array_gpu = CuArray(test_array_cpu)

CUDA.@sync reduce(*, test_array_gpu)
foldl(*, test_array_cpu)