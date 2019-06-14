% 2019 04 26  Search for parameter combination that minimizes objective


% Set paths
if ismac && isunix  % if on Mac
    addpath /Users/wu-jung/code_git/ooi_sonar/palm_nmf_matlab
    data_path = '/Users/wu-jung/code_git/ooi_sonar/sample_data/';
    save_path = '/User/wu-jung/code_git/ooi_sonar/decomp_smoothness_revisit/';
elseif isunix  % if on Linux
    addpath ~/code_git/ooi_sonar/palm_nmf_matlab
    data_path = '~/code_git/ooi_sonar/sample_data/';
    save_path = '~/code_git/ooi_sonar/decomp_smoothness_revisit';
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
sm_all = [1e3,1e6];

len = length(sm_all);
f1 = 'r';
v1 = mat2cell(3*ones(len,1),ones(len,1));
f2 = 'max_iter';
v2 = mat2cell(2000*ones(len,1),ones(len,1));
f3 = 'betaW';
v3 = mat2cell(0.1*ones(len,1),ones(len,1));
f4 = 'betaH';
v4 = mat2cell(0.1*ones(len,1),ones(len,1));
f5 = 'smoothness';
v5 = mat2cell(sm_all',ones(len,1));
params_all = struct(f1,v1,f2,v2,f3,v3,f4,v4,f5,v5);


parfor rr = 1:length(sm_all)
    fprintf('smoothness = %d\n', sm_all(rr));
    [W, H, objective, iter_times, W_init, H_init, W_steps, H_steps] = ...
        palm_nmf_detail(LL, params_all(rr));
    save_file = ...
        sprintf('%s_r%02d_betaW%2.2f_betaH%2.2f_smoothness%09.2f.mat', ...
        sname, params_all(rr).r, ...
        params_all(rr).betaW, params_all(rr).betaH, ...
        params_all(rr).smoothness);
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
