# Fast-approximation-of-Unbalanced-optimal-transport-UOT-and-Maximum-mean-discrepancy-MMD
This contribution presents significant computational accelerations to prominent schemes, which enable the comparison of measures, even with varying masses. 
Concisely, we employ nonequispaced fast Fourier transform to accelerate the radial kernel convolution in unbalanced optimal transport approximation, building on the Sinkhorn algorithm.
Accelerated schemes are presented as well for the maximum mean discrepancies involving kernels based on distances.
By employing nonequispaced fast Fourier transform, our approaches significantly reduce the arithmetic operations to compute the distances from $𝓞(n²)$ to $𝓞(n \log n)$, which enables access to large and high-dimensional data sets.

### Prerequisites

The "NFFT3.jl" package has be installed. For more details, please refer to  [https://www-user.tu-chemnitz.de](https://www-user.tu-chemnitz.de/~potts/nfft/) and https://github.com/NFFT/NFFT3.jl. 


**NOTE: The "NFFT3.jl" package should be in the same parent directory.**


### Reference

When you are using this code, please cite the paper.

**Needs to be updated**

This paper also comprehensively explains the implementation of NFFT accelerated UOT approximation and MMD.


## Directory structure

| File/Folder   | Purpose                                                                                   |
| ------------- |-------------------------------------------------------------------------------------------|   
| src           | Standard and NFFT accelerated implementations of regularized UOT and MMD  |
| Experiments | Experiment scripts-- Performance and accuracy analysis      |


