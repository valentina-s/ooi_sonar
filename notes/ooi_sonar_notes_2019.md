# Sonar Project 2019 Notes

## 2019/04/18
* Tried to figure out what parameters were used for 2017 NOLA presentation using PALM-NMF
	- It seems betaW=betaH=0.0 and eta_H=100 gives the closest results to those in the presented figures
	- The sparsity weight (lambda in the PALM-NMF code) is 0 by default; it was not changed in the 2017 test
* Need to try the cophenetic coeff and RSS approach to select the best number of components
  - see literature summary [here](https://app.nuclino.com/sonar-project/Meeting-notes/NMF-number-of-comp-review-a44403a4-c0ef-4c22-bedb-121bff378a49)
* Need to figure out a way to select eta_H
  - is it true that when eta_H is high the components are more stable? If so, can use the cophenetic coeff approach
  - how does changes in eta_H affect the residual error between the reconstruction and the data matrix? can this be used as a guide to select eta_H?
  - If we do the same comparison with permutated data matrix, do we see different slopes of RSS drop with increasing eta_H?
  - In theory we want eta_H to be large because we know the activation of patterns should be smooth across day; but when eta_H is too large, we start to lose the ability to match the daily variation in the data becomes.
    - do we see a drop and then a rise in terms of RSS when scanning across a series of eta_H?
    - there is also the question of how effective the Tihkonov regularization is in terms of forcing continuity
  - can we also quantify how similar the components (W) are as a function of eta_H?

## 2019/04/29
* Most of the work in the past week was in writing the manuscript
* On 2019/04/26 started to do parameter search over rank and smoothness in matlab.
* Today: goal is to be able to compute cophenetic coeff in matlab.
* In NIMFA: `sop`, `elop`, etc. are defined in `nimfa.utils.linalg`
* Matlab has function `cophenet` to calculate cophenetic coefficients, so have used it for checking the decomposition results in `copute_cophenetic_coeff.m`

## 2019/04/30
* Discussed with VS about how to determine the various parameters and rank
	- by changing the parameters, the optimization is actually solving different problems, so looking at how low the objective is at the end is not a meaningful measure
	- we did not understand why the residual has different slope for the shuffled data vs the original data toward higher rank, even though there is indeed a taper-off for the original data around rank 3-7 which is not observed in the shuffled data. Code:
		- `rank_search_20190426_shuffle.m`: computing smooth NMF
		- `rank_search_check_20190427_shuffle.m`: compare residual for the original and shuffled data
	- to understand this more we'll do the comparisons in the bullet points below.
* Decided that we should add PCA results as a supplementary figure with extended caption and discuss the PCA and NMF comparison in the main text.
* try comparing the residual for the original and shuffled data when using classic NMF (max iter=500):
	- `rank_search_20190430_nmf.m`: compute classic nmf for both the original and shuffled data
	- `rank_search_check_20190430_nmf_shuffle.m`: compare residual
* try adding noise and do the same comparison as the above:
	- `rank_search_20190430_nmf_Lnoise.m`: add noise, compute classic nmf
	- `rank_search_check_20190430_nmf_Lnoise_shuffle.m`: compare residual
	- it appeared that max iter=500 would not lead to convergence at higher rank, jagging it up to max iter=1000
	- even at max iter=1000 for rank>3 the decomposition for original data didn't converge
* try simulating 3 components and add noise, decompose using different rank and see what happens.
* try play with the sparsity term to see if we can get any info out of forcing sparsity


## 2019/06/14
Clean up the matlab code folder:
* ss-NMF
	- `rank_search_20190426.m`: use ss-NMF in the original `palm_nmf.m` code, sweep over different rank; results saved in `/decomp_results`; check results using `rank_search_check_20190426.m`
	- `rank_search_20190426_repeat.m`: use ss-NMF in the original `palm_nmf.m` code, sweep over different rank and make many repetitions; results saved in `/decomp_results_repeat`
	- `rank_search_check_20190426_shuffle`: use ss-NMF in the original `palm_nmf.m` code, sweep over different rank and shuffle each column; results saved in `/decomp_results`; check results using `rank_search_check_20190426_shuffle.m`; SVD was also run as a comparison case in the check.
* Classic NMF
	- `rank_search_check_20190430_nmf.m`: run classic NMF using `nnmf` in Matlab on original data and shuffled data (shuffle each column); results saved in `/decomp_results_class_nmf`; check results using `rank_search_check_20190430_nmf_shuffle.m`
	- `rank_search_check_20190430_nmf_Lnoise.m`: run classic NMF using `nnmf` in Matlab with noise-added original and shuffled data (shuffle each column); results saved in `/decomp_results_classic_nmf_Lnoise`; check results using `rank_search_check_20190430_nmf_Lnoise_shuffle.m`
	- `rank_search_20190501_nmf.m`: run classic NMF using `nnmf` in Matlab on shuffled data (shuffle each row); results saved in `/decomp_results_classic_nmf_permute_within_row`
* Normalized variance
	- `rank_fixed_20190426_repeat.m`: multiple runs of ss-NMF with a range of rank, using original data
	- `rank_fixed_20190502_repeat.m`: multiple runs of ss-NMF with a range of rank, using original data with normalized variance for each pixel (to facilitate measuring AIC)
* Cophenetic correlation coeff
	- `compute_cophenetic_coeff.m`: compute coph coeff for multiple runs with a sweep of rank
	- `compute_cophenetic_coeff_normvar.m`: compute coph coeff for multiple runs with a sweep of rank, using original data with normalized variance
	- `compute_cophenetic_coeff_normvar_modified_coeff`: testing to see if could modify the consensus matrix by using correlation between components instead of forcing 0 and 1. Unfinished.
* AIC
	- `AIC_prelim_test.m`: testing the computation
	- `AIC_classic_nmf_normalized_var.m`: compute classic NMF using original data with normalized variance, and then compute AIC
	- `AIC_ss-nmf_normalized_var.m`: compute ss-NMF using original data with normalized variance, and then compute AIC
* Smoothness factor
	- `smoothness_search_20190426.m`: sweep through a range of smoothness parameter in ss-NMF using the original data; results saved in `/decomp_results_smoothness`; check results using `smoothness_search_check_20190426.m`
	- `smoothness_search_20190611.m`: same as the above but use new function `palm_nmf_detail.m` to record initial and all steps of H and W.
* Check H and W from different runs
	- `check_H_W_multiple_runs_classic_nmf.m`
	- `check_H_W_multiple_runs_ss-nmf.m`
* Compare ss-NMF and SVD
	- `nmf_svd_cmp.m`: compare reconstruction RMS error between ss-NMF and SVD over a range of ranks
* Individual functions
	- `find_match_factor_seq.m`: function to find matching components across 2 different runs. Called by `compute_cophenetic_coeff_normvar_modified_coeff`.
	- `palm_nmf_detail.m` is new: record the initial and all updates of H and W throughout all iterations.
