function smoothness_20190621_sweep_prince_func(rank,sm,sp,beta,max_iter)

% 2019 06 21  Fix rank=3 and run a series of smoothness param
%             with a few different sparsity param
% INPUT
%   rank   rank
%   sm     smoothness
%   sp     sparsity
%   beta   betaH = betaW
%   max_iter   maximum iteration

% Set paths
addpath /home/ch153/wjl/ooi_sonar/palm_nmf_matlab
data_path = '/home/ch153/wjl/ooi_sonar/sample_data';
data_file = '20150817-20151017_MVBS_PCPcleaned.h5';
save_path = '/scratch/ch153/wjl/nmf_results/decomp_smoothness_sweep_20190621_w_lambda';

% If save_path does not exist, create it
if ~exist(save_path)
    mkdir(save_path)
end

cd(data_path)  % has to switch folder to ensure reading h5 files correctly


% Get current script name
pname = mfilename('fullpath');
[~,sname,~] = fileparts(pname);  % script name


% Load PCP-cleaned data
L = h5read(data_file,'/L');
L_sep = h5read(data_file,'/L_sep');
L_plot = h5read(data_file,'/L_plot');

depth_bin_size = h5read(data_file,'/depth_bin_size');
ping_time = h5read(data_file,'/ping_time');
ping_per_day_mvbs = h5read(data_file,'/ping_per_day_mvbs');
depth_bin_num = 37;


LL = L-min(L(:));  % make it non-negative


% Iteration to save varies by order of magnitude
iter_order = 0:4;
iter_save = repmat(1:0.2:9.8,length(iter_order),1);
for iorder = iter_order+1;
    iter_save(iorder,:) = iter_save(iorder,:)*10^iter_order(iorder);
end
iter_save = iter_save';
iter_save = unique(int32(iter_save(:)));
iter_save = iter_save(iter_save<=max_iter);
if iter_save(end)~=max_iter
    iter_save(end+1) = max_iter;
end

% Run ss-NMF
params.r = rank;
params.max_iter = max_iter;
params.betaW = beta;
params.betaH = beta;
params.smoothness = sm;
params.sparsity = sp;

fprintf('smoothness = %0.2e, sparsity = %0.2e\n',sm,sp);
[W, H, objective, iter_times, W_init, H_init, W_steps, H_steps] = ...
    palm_nmf_detail(LL, params, iter_save);
save_file = ...
    sprintf('%s_r%02d_betaW%2.2f_betaH%2.2f_smoothness%0.2e_sparsity%0.2e_maxiter%0.2e.mat', ...
            sname, params.r, ...
            params.betaW, params.betaH, ...
            params.smoothness, params.sparsity,max_iter);
save_file = fullfile(save_path, save_file);
m = matfile(save_file,'writable',true);
m.W = W;
m.H = H;
m.W_init = W_init;
m.H_init = H_init;
m.W_steps = W_steps;
m.H_steps = H_steps;
m.objective = objective;
m.iter_times = iter_times;
m.iter_save = iter_save;
m.params = params;
