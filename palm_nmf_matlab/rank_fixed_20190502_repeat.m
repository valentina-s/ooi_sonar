% 2019 05 02  Running SNMF on variance normalized data
%             with the same parameter many times


% Set paths
if ismac && isunix  % if on Mac
    addpath /Users/wu-jung/code_git/ooi_sonar/palm_nmf_matlab
    data_path = '/Users/wu-jung/code_git/ooi_sonar/sample_data/';
    save_path = '/Users/wu-jung/code_git/ooi_sonar/decomp_results_normvar_repeat/';
elseif isunix  % if on Linux
    addpath ~/code_git/ooi_sonar/palm_nmf_matlab
    data_path = '~/code_git/ooi_sonar/sample_data/';
    save_path = '~/code_git/ooi_sonar/decomp_results_normvar_repeat/';
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


% Normalize variance across pixels
LL_norm = nan(size(LL));
sigma_all = nan(size(LL,1),1);
for irow = 1:size(LL,1)
    sigma = std(LL(irow,:));
    LL_norm(irow,:) = LL(irow,:)/sigma;
    sigma_all(irow) = sigma;
end


% Run smooth NMF, sweep through different rank, repeat many times
repeat_num = 101:140;
rank_all = [1:20];

len = length(rank_all);
f1 = 'r';
v1 = mat2cell(rank_all',ones(len,1));
f2 = 'max_iter';
v2 = mat2cell(500*ones(len,1),ones(len,1));
f3 = 'betaW';
v3 = mat2cell(0.1*ones(len,1),ones(len,1));
f4 = 'betaH';
v4 = mat2cell(0.1*ones(len,1),ones(len,1));
f5 = 'smoothness';
v5 = mat2cell(100*ones(len,1),ones(len,1));
params_all = struct(f1,v1,f2,v2,f3,v3,f4,v4,f5,v5);

for nn = repeat_num
    parfor rr = 1:length(rank_all)
        fprintf('rank = %d\n', rank_all(rr));
        [W, H, objective, iter_times] = palm_nmf(LL_norm, params_all(rr));
        save_file = sprintf('%s_r%02d_betaW%2.2f_betaH%2.2f_smoothness%06.2f_rep%04d.mat', ...
                            sname, rank_all(rr), ...
                            params_all(rr).betaW, params_all(rr).betaH, ...
                            params_all(rr).smoothness, ...
                            nn);
        save_file = fullfile(save_path, save_file);
        m = matfile(save_file,'writable',true);
        m.W = W;
        m.H = H;
        m.objective = objective;
        m.iter_times = iter_times;
        m.params = params_all(rr);
    end
end