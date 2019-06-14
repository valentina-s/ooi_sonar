% 2019 04 29  hierarchical clustering and cophenetic coeff for NMF results
%             over multiple runs


if ismac && isunix  % if on Mac
    addpath /Users/wu-jung/code_git/ooi_sonar/palm_nmf_matlab/
    data_path = '/Users/wu-jung/code_git/ooi_sonar/sample_data/';
%     save_path = '/Users/wu-jung/code_git/ooi_sonar/decomp_results_repeat/';
    save_path = '/Users/wu-jung/code_git/ooi_sonar/decomp_results_normvar_repeat/';
elseif isunix  % if on Linux
    addpath ~/code_git/ooi_sonar/palm_nmf_matlab/
    data_path = '~/code_git/ooi_sonar/sample_data/';
%     save_path = '~/code_git/ooi_sonar/decomp_results_repeat/';
    save_path = '~/code_git/ooi_sonar/decomp_results_normvar_repeat/';
end


% Load PCP-cleaned data
data_file = '20150817-20151017_MVBS_PCPcleaned.h5';
cd(data_path)  % has to switch folder to ensure reading h5 files correctly
L = h5read(data_file,'/L');
L_sep = h5read(data_file,'/L_sep');
L_plot = h5read(data_file,'/L_plot');
LL = L-min(L(:));  % make it non-negative


% Normalize variance across pixels
LL_norm = nan(size(LL));
sigma_all = nan(size(LL,1),1);
for irow = 1:size(LL,1)
    sigma = std(LL(irow,:));
    LL_norm(irow,:) = LL(irow,:)/sigma;
    sigma_all(irow) = sigma;
end


% Params
smoothness = 100;
rank_all = 1:20;
fac1_rank = cell(length(rank_all),1);
fac2_rank = cell(length(rank_all),1);
AIC_rank = cell(length(rank_all),1);
AIC_rank_mean = nan(length(rank_all),1);
AIC_rank_std = nan(length(rank_all),1);
for irank = 1:length(rank_all)
    % Get filenames of results with the same param combination
    fname_pre = sprintf('%s_r%02d_betaW0.10_betaH0.10_smoothness100.00_rep*.mat', ...
        'rank_fixed_20190502_repeat', rank_all(irank));
    fname_all = dir(fullfile(save_path,fname_pre));
    
    % Concensus matrix
    AIC_ff = nan(length(fname_all),1);
    fac1 = nan(length(fname_all),1);
    fac2 = nan(length(fname_all),1);
    for ff = 1:length(fname_all)
        % Load file
        fname = fname_all(ff).name;
        D = load(fullfile(save_path, fname));
        
        % Calculate AIC
        fac1(ff) = norm(LL_norm - D.W*D.H, 'fro').^2/1.21;
        fac2(ff) = 2*rank_all(irank)*(size(LL,1)+size(LL,2));
        AIC_ff(ff) = fac1(ff) + fac2(ff);
    end
    
    fac1_rank{irank} = fac1;
    fac2_rank{irank} = fac2;
    AIC_rank{irank} = AIC_ff;
    AIC_rank_mean(irank) = mean(AIC_ff);
    AIC_rank_std(irank) = std(AIC_ff);
end