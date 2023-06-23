"""
    Maximum Mean Discrepancy (MMD) Accuracy analysis
	created: 2023, May
	author©: Rajmadan Lakshmanan
"""

using LinearAlgebra
using Distributions
using Test

@show n1=10000;n2=10000
μ= rand(Uniform(0,1),n1);μ/= sum(μ);
s1= (rand(Uniform(0,1),(n1, 2))); #s1/=sum(s1)
ν= rand(Uniform(0,1),n2);ν/= sum(ν)
s2= (rand(Uniform(0,1),n2, 2)); #s2/=sum(s2)

MMD_Gauss = MMDGauss(μ,ν,s1,s2;σ=1.0)
NFFTMMD_Gauss = NFFTMMDGauss(μ,ν,s1,s2;σ=1.0)

@show Residual_Gauss = (norm(MMD_Gauss-NFFTMMD_Gauss,2)/norm(MMD_Gauss,2)) *100

MMD_Lap = MMDLaplace(μ,ν,s1,s2;σ=1)
NFFTMMD_Lap = NFFTMMDLaplace(μ,ν,s1,s2;σ=1.0) 

@show Residual_Lap = (norm(MMD_Lap-NFFTMMD_Lap,2)/norm(MMD_Lap,2)) *100


MMD_IMQ = MMDIMQ(μ,ν,s1,s2;c=1.0)
NFFTMMD_IMQ = NFFTMMDIMQ(μ,ν,s1,s2;c=1.0) 

@show Residual_IMQ = (norm(MMD_IMQ-NFFTMMD_IMQ,2)/norm(MMD_IMQ,2)) *100


MMD_e = MMDEnergy(μ,ν,s1,s2)
NFFTMMD_e = NFFTMMDEnergy(μ,ν,s1,s2) 

@show Residual_e = (norm(MMD_e-NFFTMMD_e,2)/norm(MMD_e,2)) *100


(3.95e-12+1.31e-12+6.94e-13+1.23e-13+2.45e-12+7.28e-13+4.23e-13+1.14e-11+1.59e-11+2.21e-13)/10