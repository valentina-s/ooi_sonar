% 2019 06 17  Set rank to be higher but with sparsity constraints,
%             see if can use the results to determine rank


% Set paths
if ismac && isunix  % if on Mac
    addpath /Users/wu-jung/code_git/ooi_sonar/palm_nmf_matlab
    data_path = '/Users/wu-jung/code_git/ooi_sonar/sample_data/';
    save_path = '/User/wu-jung/code_git/ooi_sonar/decomp_smoothness_revisit_2e4';
elseif isunix  % if on Linux
    addpath ~/code_git/ooi_sonar/palm_nmf_matlab
    data_path = '~/code_git/ooi_sonar/sample_data/';
    save_path = '/media/wu-jung/internal_2tb/nmf_results/decomp_smoothness_sweep_20190621_w_lambda';
end
data_file = '20150817-20151017_MVBS_PCPcleaned.h5';

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


% Run smooth NMF, sweep through different smoothness
rank = 3;
sm_all = [1e5,2e5,5e5,1e6,2e6,5e6,1e7,2e7,5e7,1e8];
lambda_all = [10,5,2,20];
max_iter = 2e4;
beta = 0.1;

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
for ilambda = 1:length(lambda_all)

    lambda = lambda_all(ilambda);
    fprintf('sparsity = %d\n',lambda);
    
    len = length(sm_all);
    f1 = 'r';
    v1 = mat2cell(rank*ones(len,1),ones(len,1));
    f2 = 'max_iter';
    v2 = mat2cell(max_iter*ones(len,1),ones(len,1));
    f3 = 'betaW';
    v3 = mat2cell(beta*ones(len,1),ones(len,1));
    f4 = 'betaH';
    v4 = mat2cell(beta*ones(len,1),ones(len,1));
    f5 = 'smoothness';
    v5 = mat2cell(sm_all',ones(len,1));
    f6 = 'sparsity';
    v6 = mat2cell(lambda*ones(len,1),ones(len,1));
    params_all = struct(f1,v1,f2,v2,f3,v3,f4,v4,f5,v5,f6,v6);

    %parfor ilambda = 1:length(lambda_all)
    parfor ss=1:length(sm_all)
        fprintf('smoothness = %0.2e\n',sm_all(ss));
        params = params_all(ss);
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
    end

end
