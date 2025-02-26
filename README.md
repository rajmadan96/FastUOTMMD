# Unbalanced Optimal Transport and Maximum Mean Discrepancies: Interconnections and Rapid Evaluation
This contribution presents substantial computational advancements to compare measures even with varying masses.
Specifically, we utilize the nonequispaced fast Fourier transform to accelerate the radial kernel convolution in unbalanced optimal transport approximation, built upon the Sinkhorn algorithm. We also present accelerated schemes for maximum mean discrepancies involving kernels.
Our approaches reduce the arithmetic operations needed to compute distances from $𝓞(n²)$ to $𝓞(n \log n)$, opening the door to handle large and high-dimensional datasets efficiently. 
Furthermore, we establish robust connections between transportation problems, encompassing Wasserstein distance and unbalanced optimal transport, and maximum mean discrepancies.
This empowers practitioners with compelling rationale to opt for adaptable distances.

### Prerequisites

The "NFFT3.jl" package has be installed. For more details, please refer to  [https://www-user.tu-chemnitz.de](https://www-user.tu-chemnitz.de/~potts/nfft/) and https://github.com/NFFT/NFFT3.jl. 


**NOTE: The "NFFT3.jl" package should be in the same parent directory.**


### Reference

When you are using this code, please cite the paper.

<a id="1">[1]</a> Rajmadan Lakshmanan, and Alois Pichler. (2024). [Unbalanced Optimal Transport and Maximum Mean
Discrepancies: Interconnections and Rapid Evaluation](https://arxiv.org/pdf/2306.13618). 

This paper also comprehensively explains the implementation of NFFT accelerated UOT approximation and MMD.


## Directory structure

| File/Folder   | Purpose                                                                                   |
| ------------- |-------------------------------------------------------------------------------------------|   
| src           | Standard and NFFT accelerated implementations of regularized UOT and MMD  |
| Experiments | Experiment scripts-- Performance and accuracy analysis      |


