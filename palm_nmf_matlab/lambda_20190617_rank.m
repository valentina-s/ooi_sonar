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
    save_path = '/media/wu-jung/internal_2tb/nmf_results/decomp_lambda_rank_sweep_revisit_1e4';
end
data_file = '20150817-20151017_MVBS_PCPcleaned.h5';

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
rank_all = [8,6];
lambda_all = [2,5,10,20];
% The above choices came from checking results in folder
% /decomp_lambda where we can visiably see a mid-water column group
% jumped from one component to another with sparsity 5->10. Also
% sparcity=50 was too high as the components become very "empty."

sm = 1e7;
max_iter = 1e4;
beta = 0.1;

parfor rr=1:length(rank_all)

    fprintf('rank = %d\n',rank_all(rr));
    r = rank_all(rr);
    
    len = length(lambda_all);
    f1 = 'r';
    v1 = mat2cell(r*ones(len,1),ones(len,1));
    f2 = 'max_iter';
    v2 = mat2cell(max_iter*ones(len,1),ones(len,1));
    f3 = 'betaW';
    v3 = mat2cell(beta*ones(len,1),ones(len,1));
    f4 = 'betaH';
    v4 = mat2cell(beta*ones(len,1),ones(len,1));
    f5 = 'smoothness';
    v5 = mat2cell(sm*ones(len,1),ones(len,1));
    f6 = 'sparsity';
    v6 = mat2cell(lambda_all',ones(len,1));
    params_all = struct(f1,v1,f2,v2,f3,v3,f4,v4,f5,v5,f6,v6);


    for ilambda = 1:length(lambda_all)
        fprintf('lambda = %d\n', lambda_all(ilambda));
        
        params = params_all(ilambda);
        [W, H, objective, iter_times, W_init, H_init, W_steps, H_steps] = ...
            palm_nmf_detail(LL, params);
        save_file = ...
            sprintf('%s_r%02d_betaW%2.2f_betaH%2.2f_smoothness%0.2e_sparsity%0.2e.mat', ...
                    sname, params.r, ...
                    params.betaW, params.betaH, ...
                    params.smoothness, params.sparsity);
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
        m.params = params_all(rr);
    end

end
