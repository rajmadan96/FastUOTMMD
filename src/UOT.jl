"""
    Entropy Regularized Unbalanced Optimal Transport (UOT) 
	created: 2023, Febuary
	author©: Rajmadan Lakshmanan
"""

using LinearAlgebra
using LogExpFunctions


#	│	Distance function
#	╰────────────────────────────────────────────────────

function distFunction(states1::Array{Float64}, states2::Array{Float64})::Array{Float64,2}
	n1= size(states1, 1)		# works in all dimensions
	n2= size(states2, 1)
	dMatrix= Array{Float64}(undef, (n1,n2))
	for i= 1:n1, j= 1:n2		# Euclidean norm
		dMatrix[i,j]= norm(states2[j,:]- states1[i,:])
	end
	return dMatrix
end

function SinkhornUOTkl(μ::Vector{Float64}, ν::Vector{Float64}, distMatrix::Array{Float64,2}; r::Float64= 1., λ::Float64= 1.,η1::Float64= 1.,η2::Float64= 1.,tol=1e-5,max_iter=1000,verbose=false)
	βr=zeros(size(μ));
	γc= zeros(size(ν));		# guess a starting value
	count= 0; 
    distMatrix.^= r; K= exp.(-λ * distMatrix)
	while true	# Sinkhorn iteration
        βr_old = copy(βr);
        γc_old = copy(γc);

		βr= -(η1/(1+λ*η1))*(log.(vec(K*((exp.(λ*γc)).*ν))).+ 1e-20)		# vector operation
				
		γc= -(η2/(1+λ*η2))*(log.(vec(K'*((exp.(λ*βr)).*μ))).+ 1e-20)	# vector operation
		   
		count += 1;
		if count % 1 == 0
			
			err_βr=
				norm(real(βr_old)- real(βr),Inf) 
			err_γc=
				norm(real(γc_old)- real(γc),Inf) 
			if (err_βr+ err_γc) < tol ||  count > max_iter 
				break
			end
			if verbose
				println("Iteration $count, err = ",  (err_βr+ err_γc))
			end
			
		end
	end
	return (f=βr,g=γc,count=count) 
end

