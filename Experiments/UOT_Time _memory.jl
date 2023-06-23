"""
    Unbalanced Optimal Transport (UOT) -- Time and memory allocation 
	created: 2023, May
	author©: Rajmadan Lakshmanan
"""

using LinearAlgebra
using Distributions


# n and \tilde n
n1=1000000
n2=1000000

d=3                                     # dimension

μ= rand(Uniform(0,1),n1);			    # Arbitrary vector μ

s1= (rand(Uniform(0,1),(n1, d)))        # Support points x

ν= rand(Uniform(0,1),n2);               # Arbitrary vector ν

s2= (rand(Uniform(0,1),n2, d)) 	        # Support points \tilde x

r=2.0                                   # r=1 Lapalace Kernel; r=2 Gaussian Kernel;

λ=20.0 								    # Entropy regularization parameter

η=1.0;η1=η2=η                           # Marginal regularization parameter



# Standard regularied UOT algorithm

@time StdUOT= SinkhornUOTkl(μ, ν, distFunction(s1,s2); r, λ,η1,η2)

# NFFT-accelerated regularied UOT algorithm

@time NFFT_UOT = UOTSinkhornKLNFFTStabilized(μ, ν, s1,s2; r, λ,η1,η2) 

