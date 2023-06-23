"""
	NFFT Entropy Regularized Unbalanced Optimal Transport (NFFT UOT)
	created: 2023, Febuary
	author©: Rajmadan Lakshmanan
"""

using LinearAlgebra
using NFFT3
"""
	UOTSinkhornNFFT(
        μ, ν, distMatrix ; r, λ, η1,η2)
Construct an NFFT UOT Sinkhorn algorithm for computing an regularized UOT Problem.
"""

#	│	Sinkhorn-Knopp iteration algorithm , NFFT
#	╰────────────────────────────────────────────────────
function UOTSinkhornKLNFFTStabilized(μ::Vector{Float64}, ν::Vector{Float64}, S1, S2; r::Float64= 1., λ::Float64= 1.,η1::Float64= 1.,η2::Float64= 1.,tol=1e-5,max_iter=1000,verbose=false)
#	│	prepare NFFT and create a Julia plan-object
	#	╰────────────────────────────────────────────────────
    if size(μ,1) != size(S1,1) && size(ν,1) != size(S2,1) 
        error("Dimension Mismatch") 
    end
    p= 8; m= p; n= 256; flags= UInt32(0)
	eps_I= real(p/ n); eps_B= real(max(1/ 16, p/ n)); nn= 2 * n; 	
	reScale= (0.25- eps_B/ 2)/ max(1, maximum(abs.(S1)), maximum(abs.(S2)));
    if r == 1
        kernel= "laplacian_rbf"; kScale= reScale/ λ					#	exp(-λ d)
    elseif r == 2
        kernel= "gaussian";	kScale=  kScale= reScale/ sqrt(λ)		#	exp(-λ d^2)
    end
     if length(S1)==length(μ) && length(S2)==length(ν)	# one dimension
				
		plan1 = NFFT3.FASTSUM(1, length(S2), length(S1),string(kernel),real(kScale), n, p,  eps_I, eps_B, nn, m)
		plan2 = NFFT3.FASTSUM(1, size(S1,1), size(S2,1),string(kernel),real(kScale), n, p,  eps_I, eps_B, nn, m) 
		plan1.x = S2* reScale; plan1.y = S1* reScale
		plan2.x = S1* reScale; plan2.y = S2* reScale
	elseif length(S1)==2*length(μ) && length(S2)==2*length(ν)		# two dimensions
		plan1 = NFFT3.FASTSUM(2, size(S2,1), size(S1,1),kernel,real(kScale), n, p,  eps_I, eps_B, nn, m)  
		plan2 = NFFT3.FASTSUM(2, size(S1,1), size(S2,1),kernel,real(kScale), n, p,  eps_I, eps_B, nn, m) 
		plan1.x = S2* reScale; plan1.y = S1* reScale
		plan2.x = S1* reScale; plan2.y = S2* reScale
	elseif length(S1)==3*length(μ) && length(S2)==3*length(ν)		# three dimensions
		plan1 = NFFT3.FASTSUM(3, size(S2,1), size(S1,1),kernel,real(kScale), n, p,  eps_I, eps_B, nn, m)  
		plan2 = NFFT3.FASTSUM(3, size(S1,1), size(S2,1),kernel,real(kScale), n, p,  eps_I, eps_B, nn, m) 
		plan1.x = S2* reScale; plan1.y = S1* reScale
		plan2.x = S1* reScale; plan2.y = S2* reScale
    else error( "High dimension" )
	end
	#	│	iterate Sinkhorn
	#	╰────────────────────────────────────────────────────
	βr= zeros(ComplexF64, length(μ))
	γc= zeros(ComplexF64, length(ν))		# guess a starting value
	count= 0; tmpS= Vector{Float64}(undef, length(μ)); 
	while true  # Sinkhorn iteration
		βr_old = copy(βr);
        γc_old = copy(γc);
		γc= exp.(λ*γc).*ν	# rescale
		@. γc= complex(real(γc), 0.)
		tmp= reinterpret(Float64, γc)
		tmp[2:2:end].= 0.0		# force imaginary part zero
		tmp_γc= norm(γc,1)
		plan1.alpha = γc/tmp_γc;
		NFFT3.trafo(plan1)		# fast summation
		βr= -(log.(plan1.f*tmp_γc))*(η1/(1+λ*η1))
		tmp_βr= norm(exp.(λ*βr).*μ,1)
		#βr= exp.(λ*βr).*μ
		plan2.alpha= exp.(λ*βr).*μ/tmp_βr 	# vector operation
		NFFT3.trafo(plan2)
		tmpS= plan2.f
		γc= -(log.(tmpS*tmp_βr))*(η2/(1+λ*η2))			# vector operation		
		count += 1; 
		if count % 10 == 0
		err_βr=
			 norm(real(βr_old)- real(βr),Inf) 
		 err_γc=
			 norm(real(γc_old)- real(γc),Inf) 
		 if verbose
			 println("Iteration $count, err = ", (err_βr+ err_γc))
		 end
		 if (err_βr+ err_γc) < tol || count > max_iter
			 break
		 end
		end 
	 
	end 
	return (f=real(βr),g=real(γc), count= count)
end
