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
* Check H and W from different runs (try to match components by cross-correlation of H)
	- `check_H_W_multiple_runs_classic_nmf.m`
	- `check_H_W_multiple_runs_ss-nmf.m`
* Compare ss-NMF and SVD
	- `nmf_svd_cmp.m`: compare reconstruction RMS error between ss-NMF and SVD over a range of ranks
* Individual functions
	- `find_match_factor_seq.m`: function to find matching components across 2 different runs. Called by `compute_cophenetic_coeff_normvar_modified_coeff`.
	- `palm_nmf_detail.m` is new: record the initial and all updates of H and W throughout all iterations.


## 2019/06/17-19
Checking details of ss-NMF, specific questions:
* How does the objective function and H and W vary when we run many many iterations? How does H and W vary throughout the iterations? --> check for convergence issues
* Can we try to push sparsity but decompose at a higher rank to see if we can discover dimensionality that way? i.e., only a few components will be nonzero.
* How does making betaW=betaH=0 change the iterative updates of objective function and the resulting W and H components.

Code added to investigate the above:
* beta
	- `beta_20190616.m`: run ss-NMF with the following params:
		```~matlab
		rank = 3;
		smoothness: [1e2, 1e4, 1e6, 1e8];
		max_iter = 1e4;
		betaW = 0;
		betaH = 0;
		```
		- `/decomp_beta_revisit_1e4` has results from 1e4 iterations
		- `/decomp_beta` has results from 2e3 iterations
	- `check_beta.m`: check results from the above by plotting W and H at various iterations and overall changes of the objective function
* lambda (sparsity)
	- `lambda_20190616.m`: run ss-NMF with the following params:
		```~matlab
		smoothness = 1e6;
		max_iter = 1e3;
		betaW = betaH = 0.1
		sparsity = [1,2,5,10,20,50];
		```
		results saved in `/decomp_lambda`
	- `lambda_20190617_rank.m`: run ss-NMF with high rank (>3, e.g., 10, 8, ...) but push for sparsity to see if can find out the approximate rank of the data matrix (i.e., true dimension will be the number of nonzero Ws).
		```~matlab
		smoothness = 1e7;
		max_iter = 1e4;
		betaW = betaH = 0.1;
		sparsity = [2,5,10,20];
		rank = [10, 8, 6];   % rank 6 and 8 not completed yet
		```
		results saved in `/decomp_lambda_rank_sweep_revisit_1e4`:
	- `lambda_20190620_rank_large_step.m`: run ss-NMF with rank=10 and push for sparsity; this is to check if the nonzero Ws from `lambda_20190617_rank` are stable when we run many many iterations. Parameters:
		```~matlab`
		smoothness = 1e7;
		max_iter = 5e4;
		betaW = betaH = 0.1;
		sparsity = [2,5,10,20];
		rank = 10;
		```
	- `check_lambda.m`: check results from `lambda_20190616.m`
	- `check_lambda_rank_sweep_revisit_1e4.m`: check results from `lambda_20190617_rank`

## 2019/06/20
* Modified `palm_nmf_detail.m` so that can save only a subset of iterations by passing an varargin. The code is backward compatible, so with no additional varargin all iterations are saved. This was because the memory problem for ss-NMF decomposition at higher rank is due to the large matrix that needs to be kept around when all iterations are saved (--> note this is because the results are only saved at once when the whole decomposition is finished. Consider changing it so that each step is saved incrementally --> can we do this for MAT files??).
* Added `subsample_saved_results.m` to subsample previously saved results to reduce space requirements, especially since at large iterations the objective varies only very slowly. The subsampling takes equal numbers of points for each order of magnitude.
	- reduced file size in folders:
		- `/decomp_beta`
		- `/decomp_lambda`
		- `/decomp_lambda_rank_sweep`
		- `/decomp_lambda_rank_sweep_revisit_1e4` --> accidentally didn't save the first 1:9 steps in the subsampled data due to programming errors... on well.
		- `/decomp_smoothness_revisit_2e4`
* Updated W/H/objective checking code to work with subsampled data:
	- `check_lambda_rank_sweep_revisit_1e4.m`
	- `check_smoothness.m`
	- `check_beta.m`
* Added code `smoothness_20190621_sweep.m` to run ss-NMF using a series of smoothness parameters. Wanted to see if it is possible to come up with a quantitative way to select smoothness? for example, check the resulting H vectors and see if there is a jump of some sort while using some kind of smoothness criteria.
* Trying to run parallel job *embarrassingly* through array job on NYU cluster.
* Temp summary from various param searches above:
	- The smoothness effects start to show with smoothness param set to 1e6 and the H curves become quite smooth at 1e8. But even then the structures from the initial non-smoothed curves are preserved, which is a good thing. Now doing a sweep of smoothness param to determine which values to take. Rank = 3 for these runs.
	- Pushing for sparsity indeed produced the effect of zeroing out most Ws and leave only a few nonzero ones:
		- with rank=3, lambda=100 one of the 3 components is already zero. lambda=10 seems a little better.
			- These are results in folder `/decomp_lambda` and was run for only 1000 iterations. At higher iterations the sparsity effects may be stronger. These results give a good guide on selecting maybe lambda<100 for test with higher rank.
			- We can visibly see a mid-water column group jumped from one component to another with lambda 5->10. Also lambda=50 was too high as the components become very *empty*.
		- with `rank=10, lambda=[2,5,10,20], max_iter=1e4`, we can see in the folder `/decomp_lambda_rank_sweep_revisit_1e4`:
			- `lambda=2`: almost all Ws have some nonzero pixels but it is visibly more sparse than when lambda=0.
			- `lambda=5`: ~5 Ws have strong nonzero pixels and ~3 have weaker nonzero pixels.
			- `lambda=10`: only 3 Ws have nonzero pixels
			- `lambda=20`: only 3 Ws have nonzero pixels, the nonzero Ws are similar to those left in sparsity=10.
			- I suspect setting lambda somewhere around 10-20 may be the way to go for final analysis.
			- the objective functions do not have oscillation as those seen with lambda=0, which is a good thing.
TODO:
- need to subsample the `decomp_lambda` folder
- check what's in the `decomp_lambda_rank_sweep` folder


## 2019/06/23-24
- Was able to send array jobs to prince hpc cluster
	- matlab function called is `ssNMF_runner.m`
	- scripts used to send jobs are called `run_ssNMF_array_20190622_rank**.sh`
	- script `generate_ssNMF_param_file.m` is used to generate the parameter files for running array jobs
	- folder `/scratch/ch153/wjl/ssNMF_runner_param_files_20190623` contains parameter files that were used to run array jobs for rank 3-10, smoothness order 5,6,7,8, and sparsity 1,2,5,10,20,50
	- results are saved in `/ssNMF_sweep_sm_sp_20190623`
- add code `check_sm_sp_sweep.m` to plot H and W at specific iterations, so that we can check how varying smoothness and sparsity parameters would change the results. Specifically we want to see
	1. How many W components get pushed to 0 --> may give us some hint about the dimensionality of the data, also will help us select the appropriate sparsity parameter
	2. Can we use these results to come up with a measure/metric for selecting the smoothness parameter?

## 2019/06/26
- Send new array jobs to run for: rank 3-10, smoothness order 1,2,3,4, sparsity 1,2,5,10,20,50. (Only the smoothness order was different from the 20190623 runs)
- results are saved in `/ssNMF_sweep_sm_sp_20190626`

## 2019/08/02
(Revamping this work after a whole month of travel in July!)

- rsync smoothness/sparsity results from NYU HPC to MURI_4TB hard drive
	```~shell
	rsync -avzhe ssh ch153@prince.hpc.nyu.edu:/scratch/ch153/wjl/nmf_results/ssNMF_sweep_sm_sp_20190623
	rsync -avzhe ssh ch153@prince.hpc.nyu.edu:/scratch/ch153/wjl/nmf_results/ssNMF_sweep_sm_sp_20190626
	```
	confirmed that all files are downloaded
- use `check_sm_sp_sweep.m` to plot all H and W and objective curve for all smoothness and sparsity combinations, results saved in folder `ssNMF_sweep_sm_sp_check` under MURI_4TB/nmf_results
- set up NYU HPC to run sparsity=0 case
- how to decide what smoothness parameter to use?
	- `select_smoothness.m` calculates smoothness score (= lag-1 autocorrelation coefficient) of the *normalized* time series of H.
	- the tapering off of smoothness score and reconstruction error is ~sm=5e7, depending on the sparsity parameter
	- however when checking the actual W results, it seems that sm=5e7 produces too sparse a result and sm=5e6 or sm=1e7 is more reasonable.
- how to decide what sparsity parameter to use?
	- ended up plotting the location of where the minimum objective is as a function of sm and sp combinations
	- `select_sparsity.m`
	- it seems that there is a zoom where ssNMF would have oscillatory behavior and not monotonically decreasing with increasing iterations --> not sure why it is happening as PALM is supposed to be monotonic?


## 2019/09/13
(Revamping this work after another month of craziness in August!)

- downloaded results from 20190802 to `MURI_4TB/nmf_results/ssNMF_sweep_sm_sp_20190802` (sparsity=0 case)
- set up another batch of repeated runs with the same parameter combination to check for stability when both sparsity and smoothness are constrained.
	- the runs call `ssNMF_runner_repeat` (modified from `ssNMF_runner` by adding _rep%03d_ at the end of the filename)
	- rank=3, smoothness=5e6, sparsity=5 (`generate_ssNMF_param_file_20190913_r03_sm5e6_sp5e0.m`)
	- rank=4, smoothness=5e6, sparsity=10 (`generate_ssNMF_param_file_20190913_r04_sm5e6_sp1e1.m`)
	- rank=5, smoothness=1e7, sparsity=10 (`generate_ssNMF_param_file_20190913_r05_sm1e7_sp1e1.m`)
	- rank=10, smoothness=5e6, sparsity=10 (`generate_ssNMF_param_file_20190913_r10_sm5e6_sp1e1.m`)
	- the results are stored under `MURI_4TB/nmf_results/ssNMF_sweep_sm_sp_20190913`


## 2019/09/15
- set up repeated runs to find rank, using the same smoothness and sparsity combinations
	- params sent: rank=2-10, smoothness=5e6, sparsity=5
	- to generate parameter files: `generate_ssNMF_param_file_20190913_r03_sm5e6_sp5e0.m`
	- the runs call `ssNMF_runner_normvar_repeat` (modified from `ssNMF_runner` by adding _rep%03d_ at the end of the filename and normalized the data matrix within each pixel before the decomposition)
	- use sbatch to send scripts `run_ssNMF_array_20190915_rankXX.sh` (one script per rank)
- temporarily renamed `ssNMF_sweep_sm_sp_20190913` to `ssNMF_sweep_sm_sp_20190913_old` because forgot to change the path in `ssNMF_runner_normvar_repeat`


## 2019/10/15
- started revising tensor notebook for OCEANS'19 conference paper
- cannot reproduce the same non-negative parafac results as in 2018 Victoria ASA using `tensorly==0.4.4`
- looked back to lab desktop and can reproduce results as in 2018 Victoria ASA
	- the environment is `test_nmf`, with `tensorly==0.4.2`
	- it appeared that the results before might not have converged
- now using `tensorly==0.4.4` and set `tol=1e-10` to force decomposition to run very long
	- initialized with svd
	- it seems that the results are more interesting too:
		- components 0 and 1 have similar frequency signature, but different daily echogram patterns and activation patterns
		- components 1 and 2 have similar echogram pattern, but opposite frequency trends and activation patterns


## 2019/10/20-24
- In tensorly, the function `non_negative_parafac` now is actually just a wrapper over `parafac` but lacks a few functions. Therefore changed all calls to use `parafac` but setting `non_negative=True` in the intpus.
- Found out that by setting optional parameter `normalization_factor=True`, the NTF results will converge regardless of the initial condition (`init='svd'` or `init='random'`)
- Now using `svd` to ensure the order of components are the same in each run.
- Quick notebook for checking init conditions: `2019-10-21_compare-NTF-init.ipynb`
- 2019 OCEANS analysis notebook: `2019-10-22_oceans-conf-paper-analysis.ipynb`
- 2019 OCEANS plotting notebook: `2019-10-23_oceans-conf-paper-plot.ipynb`


# 2019/10/25
- Found out that the SlurmArray jobs were all started using the same initial condition, so the repetitions are identical... have to set it up right and run again.
- Turned out that on cluster, when each Matlab job is spun up, (from [stackoverflow](https://stackoverflow.com/questions/49097078/randi-generating-the-same-result-on-a-cluster))

  > the nodes you're getting the same values because the cluster nodes each act on a new instance on MATLAB, which uses rng('default') as the initialisation for the random number generator.

- The easiest way to solve this is to set `rng(row_num)` since we just need each element to be different across separate runs, so that off-by-one increament of the rng is good enough to achieve that.


# 2019/10/26
- For some reason it seems that the `rank=3, sm=5e6, sp=5` results from 20190623 batch of ssNMF calculation is a special case, in that there are 3 components left even at step 20000. For almost all the runs with the same parameter combination from 20191025, the 3rd component vanishes after ~step 2000.
- Confirmed that the difference induced by the different initial conditions from 20190623 run and 20191025 runs are real: can reproduce the same results by running locally using the 20190623 initial W and H conditions. This is in folder `ssNMF_repeat_20191026_runaway` and files:
 	- `20190623_ssNMF_runner_normvar_repeat_r05_betaW0.10_betaH0.10_smoothness5.00e+06_sparsity5.00e+00_maxiter1.20e+04_rep006.mat`

  - `20191025_ssNMF_runner_normvar_repeat_r05_betaW0.10_betaH0.10_smoothness5.00e+06_sparsity5.00e+00_maxiter2.00e+03_rep006.mat`

- Issued another round of runs for `sparsity=1` in folder `ssNMF_repeat_20191026` for rank=3-10, sm=5e6
  - The results seem less stable compare to setting sparsity=5: approximately 4 components would be nonzero at step 20000, with 1 component almost 0. But the 4 components are different from run to run. In the case when sparsity=5, in general only 2 components will remain nonzero, and the 2 are very similar.


# 2019/10/27
- Try to see if controlling betaH and betaW can prevent the runaway behavior.
  - Use the same parameter as the 20190623 file `ssNMF_runner_normvar_repeat_r05_betaW0.10_betaH0.10_smoothness5.00e+06_sparsity5.00e+00_maxiter1.20e+04_rep006.mat` but set betaH and betaW to 1.
	- The results will be in folder `ssNMF_repeat_20191026_runaway`, file: `20190623_ssNMF_runner_normvar_repeat_r05_betaW1.00_betaH1.00_smoothness5.00e+06_sparsity5.00e+00_maxiter1.20e+04_rep006.mat` (note betaH=**1.00** and betaW=**1.00**)


# 2019/10/28
- Results from varying betaH and betaW: it seems that by pushing these 2 params much larger, the decomposition will be dominated by the related terms and then won't be able	



# Summary of folders
- `ssNMF_repeat_20190623`: rank=2-10, sm=[1,2,5]x[1e5,1e6,1e7,1e8], sp=[1,2,5,10,20,50], max_iter=2e4, betaW=betaH=0.1
- `ssNMF_repeat_20190626`: rank=2-10, sm=[1,2,5]x[1e1,1e2,1e3,1e4], sp=[1,2,5,10,20,50], max_iter=2e4, betaW=betaH=0.1
- `ssNMF_repeat_20191025`: rank=2-10, sm=5e6, sp=5, max_iter=2e4, betaW=betaH=0.1
- `ssNMF_repeat_20191026`: rank=2-10, sm=5e6, sp=1, max_iter=2e4, betaW=betaH=0.1
- `ssNMF_repeat_20191027`: rank=2-10, sm=5e6, sp=2, max_iter=2e4, betaW=betaH=0.1
- `ssNMF_betaHW_20191027`: rank=2-10, sm=5e6, sp=2, max_iter=2e3, betaW=betaH=[0, 0.1, 1, 10, 100, 1000, 10000]
