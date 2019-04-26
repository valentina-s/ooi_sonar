% 2017 11 27  Test PALM-NMF algorithm on MVBS data

% clear

addpath /Users/wu-jung/code_git/ooi_sonar/

data_path = '/Users/wu-jung/code_git/ooi_sonar/sample_data/';
% save_path = '/media/wu-jung/wjlee_apl_2/ooi_zplsc_palm_nmf_2018/';
data_file = '20150817-20151017_MVBS_PCPcleaned.h5';

cd(data_path)  % has to switch folder to ensure reading h5 files correctly
%h5disp(data_file)


% Load PCP-cleaned data
L = h5read(data_file,'/L');
L_sep = h5read(data_file,'/L_sep');
L_plot = h5read(data_file,'/L_plot');
% S = h5read(data_file,'/S');
% S_sep = h5read(data_file,'/S_sep');
% S_plot = h5read(data_file,'/S_plot');
depth_bin_size = h5read(data_file,'/depth_bin_size');
ping_time = h5read(data_file,'/ping_time');
ping_per_day_mvbs = h5read(data_file,'/ping_per_day_mvbs');
depth_bin_num = 37;


LL = L-min(L(:));  % make it non-negative


% Run PALM-NMF
params.r = 3;
params.max_iter = 500;
params.betaW = 1;
params.betaH = 1;
params.smoothness = 1000;
[W, H, objective, iter_times] = palm_nmf(LL, params);


