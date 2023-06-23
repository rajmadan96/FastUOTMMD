"""
	Maximum Mean Discrepancies (MMDs)
	created: 2023, Febuary
	author©: Rajmadan Lakshmanan
"""

using LinearAlgebra

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
#	│	Prominent Kernel functions
#	╰────────────────────────────────────────────────────

function GaussKernel(x::Array{Float64},y::Array{Float64}; σ=1.)	#	Gauss kernel
	return exp.(-((distFunction(x,y)).^2/ σ^2 ))
end

function LaplaceKernel(x::Array{Float64},y::Array{Float64}; σ=1.)	#	Laplace kernel
	return exp.(-(distFunction(x,y)/ σ ))
end

function IMQKernel(x::Array{Float64},y::Array{Float64}; c=1.)	#	Inverse multiquadric kernel
	return ((distFunction(x,y)).^2 .+ c^2 ).^(-1/2)
end

function EnergyKernel(x::Array{Float64},y::Array{Float64})	#	Energy kernel
	return(-distFunction(x,y))
end

#	│	Standard MMD 
#	╰────────────────────────────────────────────────────
function MMDGauss(α::Vector{Float64},β::Vector{Float64}, x::Array{Float64},y::Array{Float64}; σ=1.)#  MMD Gauss 
    k1= GaussKernel(x,x; σ)
    k2= GaussKernel(x,y; σ)
    k3= GaussKernel(y,y; σ)
    return dot(α,k1*α) -2dot(α,k2*β)+dot(β,k3*β)
end

function MMDLaplace(α::Vector{Float64},β::Vector{Float64}, x::Array{Float64},y::Array{Float64}; σ=1.)#    MMD Laplace 
    k1= LaplaceKernel(x,x; σ)
    k2= LaplaceKernel(x,y; σ)
    k3= LaplaceKernel(y,y; σ)
	return dot(α,k1*α) -2dot(α,k2*β)+dot(β,k3*β)
end

function MMDIMQ(α::Vector{Float64},β::Vector{Float64}, x::Array{Float64},y::Array{Float64};c=1.)#     MMD Inverse multiquadric 
    k1= IMQKernel(x,x; c)
    k2= IMQKernel(x,y; c)
    k3= IMQKernel(y,y; c)
    #@show k2*β
	return dot(α,k1*α)-2dot(α,k2*β)+dot(β,k3*β)
end

function MMDEnergy(α::Vector{Float64},β::Vector{Float64}, x::Array{Float64},y::Array{Float64})#	MMD Energy (only probability measures)
    summ_α = sum(α); summ_β = sum(β)
    if summ_α != 1
        α= α./summ_α
        @info "Given vector α is not a probability vector. Thus, normalized." 
    end
    if summ_β != 1
        β= β./summ_β
        @info "Given vector β is not a probability vector. Thus, normalized." 
    end
    k1= EnergyKernel(x,x)
    k2= EnergyKernel(x,y)
    k3= EnergyKernel(y,y)
    	return dot(α,k1*α) -2dot(α,k2*β)+dot(β,k3*β)
end
