"""
    Unbalanced Optimal Transport (UOT) -- Accuracy analysis
	created: 2023, May
	author©: Rajmadan Lakshmanan
"""

using LinearAlgebra
using Distributions
using Images
using FileIO
using DataFrames


img_path = "/HOME1/users/personal/lraj/Nextcloud/Alois & Raj/Sinkhorn with Maringals/CODE/DOTMARK/WhiteNoise/picture64_1001.png"
img_0 = load(img_path)
#img_0 = imresize(img_0,(32,32))
img_0 = imresize(img_0,(64,64))
#img_0 = imresize(img_0,(128,128))
(n1,n2)= size(img_0)
s1= Array{Float64}(undef, (n1*n2,2)) # Support points x
for i=1: n1
	for j=1:n2
		s1[(i-1)*n2+j,1]= i/n1
		s1[(i-1)*n2+j,2]= j/n2
end;end

imgg_0= Gray.(img_0)
mat_0 = convert(Array{Float64}, imgg_0)
μ = vec(mat_0)  # Arbitrary vector μ

img_path_2 = "/HOME1/users/personal/lraj/Nextcloud/Alois & Raj/Sinkhorn with Maringals/CODE/DOTMARK/WhiteNoise/picture64_1002.png"
img_2 = load(img_path_2)
#img_2 = imresize(img_2,(32,32))
img_2 = imresize(img_2,(64,64))
#img_2 = imresize(img_2,(128,128))
(n1,n2)= size(img_2)
s2= Array{Float64}(undef, (n1*n2,2)) # Support points \tilde x
for i=1: n1
	for j=1:n2
		s2[(i-1)*n2+j,1]= i/n1
		s2[(i-1)*n2+j,2]= j/n2
end;end

imgg_2= Gray.(img_2)
mat_2 = convert(Array{Float64}, imgg_2)
ν = vec(mat_2) # Arbitrary vector ν


r=2.0                                   # r=1 Lapalace Kernel; r=2 Gaussian Kernel;

λ=20.0 								    # Entropy regularization parameter

η=1.0;η1=η2=η                           # Marginal regularization parameter



# Standard regularied UOT algorithm

@time StdUOT= SinkhornUOTkl(μ, ν, distFunction(s1,s2); r, λ,η1,η2)

# NFFT-accelerated regularied UOT algorithm

@time NFFT_UOT = UOTSinkhornKLNFFTStabilized(μ, ν, s1,s2; r, λ,η1,η2) 

β= StdUOT.β; γ= StdUOT.γ

β_NFFT =NFFT_UOT.f; γ_NFFT =NFFT_UOT.g

# Relative error in %
@show ((norm(β .- β_NFFT,2)/norm(β,2)) +  (norm(γ .- γ_NFFT,2)/norm(γ,2)) )*100