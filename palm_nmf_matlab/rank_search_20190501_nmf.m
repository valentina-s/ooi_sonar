% 2019 04 30  Run classic NMF over different rank


% Set paths
if ismac && isunix  % if on Mac
    addpath /Users/wu-jung/code_git/ooi_sonar/palm_nmf_matlab
    data_path = '/Users/wu-jung/code_git/ooi_sonar/sample_data/';
    save_path = '/Users/wu-jung/code_git/ooi_sonar/decomp_results_classic_nmf_permute_within_row/';
elseif isunix  % if on Linux
    addpath ~/code_git/ooi_sonar/palm_nmf_matlab
    data_path = '~/code_git/ooi_sonar/sample_data/';
    save_path = '~/code_git/ooi_sonar/decomp_results_classic_nmf_permute_within_row/';
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


% Shuffle observation
rng(168)
rand_seq = nan(size(LL));
LL_perm = nan(size(LL));
for ii=1:size(LL,1)  % shuffle within row
    rand_seq(ii,:) = randperm(size(LL,2));
    LL_perm(ii,:) = LL(ii,rand_seq(ii,:));
end


% Run smooth NMF, sweep through different rank
rank_all = 1:30;
opt = statset('MaxIter',2000,'Display','final','TolFun',1e-6);

parfor rr = 1:length(rank_all)
    fprintf('rank = %d\n', rank_all(rr));
    
    % Original data
    [W, H, D] = nnmf(LL,rank_all(rr),'option',opt,'algorithm','mult');
    save_file = sprintf('%s_r%02d.mat', ...
        sname, rank_all(rr));
    save_file = fullfile(save_path, save_file);
    m = matfile(save_file,'writable',true);
    m.W = W;
    m.H = H;
    m.D = D;

    % Shuffled data
    [W, H, D] = nnmf(LL_perm,rank_all(rr),'option',opt,'algorithm','mult');
    save_file = sprintf('%s_r%02d_shuffle.mat', ...
        sname, rank_all(rr));
    save_file = fullfile(save_path, save_file);
    m = matfile(save_file,'writable',true);
    m.W = W;
    m.H = H;
    m.D = D;
    m.rand_seq = rand_seq;
end
