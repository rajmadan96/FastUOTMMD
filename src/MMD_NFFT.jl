"""
	NFFT Maximum Mean Discrepancies (NFFT MMDs)
	created: 2023, Febuary
	author©: Rajmadan Lakshmanan
"""

using NFFT3
using LinearAlgebra


function NFFTkernel(μ::Vector{Float64}, ν::Vector{Float64}, S1::Array{Float64}, S2::Array{Float64}, kernel::String; c::Float64=1.)
    #	│	prepare NFFT and create a Julia plan-object
        #	╰────────────────────────────────────────────────────
        if size(μ,1) != size(S1,1) && size(ν,1) != size(S2,1) 
            error("Dimension Mismatch")
        end
        p= Int32(8); m= p; n= Int32(256); flags= UInt32(0)
	    eps_I= real(p/ n); eps_B= real(max(1/ 16, p/ n)); nn= 2 * n; 	#	exp(-d^2/σ^2)
	    reScale= (0.25- eps_B/ 2)/ max(1, maximum(abs.(S1)), maximum(abs.(S2)));
        if kernel == "gaussian"
            kScale= reScale*c #	exp(-d^2/c^2)
        elseif kernel == "laplacian_rbf"
            kScale= reScale*c #	exp(-d/c)
        elseif kernel == "inverse_multiquadric"
            kScale= reScale*c #	1/sqrt(d^2+c^2)
        elseif kernel == "absx"
            kScale= 1/reScale  #	|d|
        end
         if length(S1)==length(μ) && length(S2)==length(ν)	# one dimension
            plan1 = NFFT3.FASTSUM(1, length(S2), length(S1),kernel,real(kScale), n, p,  eps_I, eps_B, nn, m)
            plan1.x = S2* reScale; plan1.y = S1* reScale
        elseif length(S1)==2*length(μ) && length(S2)==2*length(ν)		# two dimensions
            plan1 = NFFT3.FASTSUM(2, size(S2,1), size(S1,1),kernel,real(kScale), n, p,  eps_I, eps_B, nn, m)
            plan1.x = S2* reScale; plan1.y = S1* reScale
        elseif length(S1)==3*length(μ) && length(S2)==3*length(ν)		# three dimensions
            plan1 = NFFT3.FASTSUM(3, size(S2,1), size(S1,1),kernel,real(kScale), n, p,  eps_I, eps_B, nn, m)#,Int32(1),Int32(1),flags
            plan1.x = S2* reScale; plan1.y = S1* reScale
        else error( "High dimension" )
        end
        ν=Vector{ComplexF64}(ν) # prerequisite
        @. ν= complex(real(ν), 0.)
        tmp= reinterpret(Float64, ν)
        tmp[2:2:end].= 0.0		# force imaginary part zero
        tmp= norm(ν,2) # for normalization
        tmp_ν= ν/tmp
        plan1.alpha = tmp_ν; NFFT3.trafo(plan1);
        if kernel == "inverse_multiquadric"
        result =real(plan1.f)*(reScale)*(tmp)
        elseif kernel == "absx" 
            result =real(plan1.f)*(1/reScale)*(tmp) # absx no need normalization because the vector here is a probability vector
        else
            result =real(plan1.f)*(tmp)
        end;
        return result
    end
#	│	NFFT MMD (Gauss kernel); 1, 2 and 3 dimension 
#	╰────────────────────────────────────────────────────
function NFFTMMDGauss(μ::Vector{Float64}, ν::Vector{Float64}, x::Array{Float64}, y::Array{Float64}; σ::Float64=1.)
	#	│	prepare three NFFT fast summation
	#	╰────────────────────────────────────────────────────
    
    s1= NFFTkernel(μ,μ,x,x,"gaussian"; c=σ)
    s2= NFFTkernel(μ,ν,x,y,"gaussian"; c=σ)
    s3= NFFTkernel(ν,ν,y,y,"gaussian";c=σ)
	return dot(μ,real(s1)) -2dot(μ,real(s2))+dot(ν,real(s3))
end

#	│	NFFT MMD (Laplace kernel); 1, 2 and 3 dimension
#	╰────────────────────────────────────────────────────
function NFFTMMDLaplace(μ::Vector{Float64}, ν::Vector{Float64}, x::Array{Float64}, y::Array{Float64}; σ::Float64=1.)
	#	│	prepare three NFFT fast summation
	#	╰────────────────────────────────────────────────────
    #c=σ
    s1= NFFTkernel(μ,μ,x,x,"laplacian_rbf"; c=σ)
    s2= NFFTkernel(μ,ν,x,y,"laplacian_rbf"; c=σ)
    s3= NFFTkernel(ν,ν,y,y,"laplacian_rbf";c=σ)
	return dot(μ,real(s1)) -2dot(μ,real(s2))+dot(ν,real(s3))
end

function NFFTMMDIMQ(μ::Vector{Float64}, ν::Vector{Float64}, x::Array{Float64}, y::Array{Float64}; c::Float64=1.)
	#	│	prepare three NFFT fast summation
	#	╰────────────────────────────────────────────────────
  
    s1= NFFTkernel(μ,μ,x,x,"inverse_multiquadric"; c)
    s2= NFFTkernel(μ,ν,x,y,"inverse_multiquadric"; c)
    s3= NFFTkernel(ν,ν,y,y,"inverse_multiquadric";c)
	return dot(μ,real(s1)) -2dot(μ,real(s2))+dot(ν,real(s3))
end

#	│	NFFT MMD (Energy kernel); 1, 2 and 3 dimension
#	╰────────────────────────────────────────────────────
function NFFTMMDEnergy(μ::Vector{Float64}, ν::Vector{Float64}, x::Array{Float64}, y::Array{Float64})
    summ_μ = sum(μ); summ_ν = sum(ν)
    if summ_μ != 1
        μ= μ./summ_μ
        @info "Given vector μ is not a probability vector. Thus, normalized." 
    end
    if summ_ν != 1
        ν= ν./summ_ν
        @info "Given vector ν is not a probability vector. Thus, normalized." 
    end
    #	│	prepare three NFFT fast summation
	#	╰────────────────────────────────────────────────────
    s1= NFFTkernel(μ,μ,x,x,"absx")
    s2= NFFTkernel(μ,ν,x,y,"absx")
    s3= NFFTkernel(ν,ν,y,y,"absx")
	return -dot(μ,real(s1)) +2dot(μ,real(s2))-dot(ν,real(s3))
end
