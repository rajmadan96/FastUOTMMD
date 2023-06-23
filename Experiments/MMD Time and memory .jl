"""
    Maximum Mean Discrepancy (MMD) Time and memory test
	created: 2023, May
	author©: Rajmadan Lakshmanan
"""

using LinearAlgebra
using Distributions
using Test

@show n1=1000;n2=1000
μ= rand(Uniform(0,1),n1);#μ/= sum(μ);
s1= (rand(Uniform(0,1),(n1, 3)));
ν= rand(Uniform(0,1),n2);#ν/= sum(ν)
s2= (rand(Uniform(0,1),n2, 3));

@time "Time and memory Gaussian Kernel" MMDGauss(μ,ν,s1,s2;σ=1.0)
@time "Time and memory NFFT Gaussian Kernel" NFFTMMDGauss(μ,ν,s1,s2;σ=1.0)


@time "Time and memory Laplace Kernel" MMDLaplace(μ,ν,s1,s2;σ=1)
@time "Time and memory NFFT Laplace Kernel" NFFTMMDLaplace(μ,ν,s1,s2;σ=1.0) 


@time "Time and memory IMQ Kernel" MMDIMQ(μ,ν,s1,s2;c=5.0)
@time "Time and memory NFFT IMQ Kernel" NFFTMMDIMQ(μ,ν,s1,s2;c=1.0) 


@time "Time and memory Energy Kernel" MMDEnergy(μ,ν,s1,s2)
@time "Time and memory NFFT Energy Kernel" NFFTMMDEnergy(μ,ν,s1,s2) 


