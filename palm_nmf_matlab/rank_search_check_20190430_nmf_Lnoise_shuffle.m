% 2019 04 26  Check param search results


if ismac && isunix  % if on Mac
    addpath /Users/wu-jung/code_git/ooi_sonar/palm_nmf_matlab/
    data_path = '/Users/wu-jung/code_git/ooi_sonar/sample_data/';
    save_path = '/Users/wu-jung/code_git/ooi_sonar/decomp_results_classic_nmf_Lnoise/';
elseif isunix  % if on Linux
    addpath ~/code_git/ooi_sonar/palm_nmf_matlab/
    data_path = '~/code_git/ooi_sonar/sample_data/';
    save_path = '~/code_git/ooi_sonar/decomp_results_classic_nmf_Lnoise/';
end

% Load PCP-cleaned data
data_file = '20150817-20151017_MVBS_PCPcleaned.h5';
cd(data_path)  % has to switch folder to ensure reading h5 files correctly
L = h5read(data_file,'/L');
L_sep = h5read(data_file,'/L_sep');
L_plot = h5read(data_file,'/L_plot');
LL = L-min(L(:));  % make it non-negative

% Sweep params
smoothness = 10;
rank_all = 1:30;

RE_original = nan(length(rank_all),1);
RE_shuffle = nan(length(rank_all),1);
for rr = 1:length(rank_all)
    % original data
    fname = sprintf('%s_r%02d.mat', ...
        'rank_search_20190430_nmf_Lnoise', rank_all(rr));
    D_original = load(fullfile(save_path, fname));
    RE_original(rr) = norm(LL - D_original.W*D_original.H, 'fro');
    
    % shuffled data
    fname = sprintf('%s_r%02d_shuffle.mat', ...
        'rank_search_20190430_nmf_Lnoise', rank_all(rr));
    D_shuffle = load(fullfile(save_path, fname));
    LL_shuffle = nan(size(LL));
    for icol = 1:size(LL_shuffle,2)
        LL_shuffle(:,icol) = LL(D_shuffle.rand_seq(:,icol),icol);
    end
    RE_shuffle(rr) = norm(LL_shuffle - D_shuffle.W*D_shuffle.H, 'fro');
end