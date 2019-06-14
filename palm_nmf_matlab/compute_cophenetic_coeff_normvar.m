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


% Params
smoothness = 100;
rank_all = 1:20;
obs_len = size(LL,2);
c_rank = [];
c_rank_weighted = [];
for rank = rank_all
    
    % Get filenames of results with the same param combination
    fname_pre = sprintf('%s_r%02d_betaW0.10_betaH0.10_smoothness100.00_rep*.mat', ...
        'rank_fixed_20190502_repeat', rank);
    fname_all = dir(fullfile(save_path,fname_pre));
    fname_all = fname_all(1:130);
    
    % Concensus matrix
    C = nan(obs_len, obs_len, length(fname_all));
    RE_file = nan(length(fname_all),1);
    for ff = 1:length(fname_all)
        % Load file
        fname = fname_all(ff).name;
        load(fullfile(save_path, fname));
        
        % Get RE
        RE_file(ff) = norm(LL_norm - W*H, 'fro');
        
        % Connectivity matrix
        [~, Hm_idx] = max(H, [], 1);  % membership assignment
        HmHm = repmat(Hm_idx, obs_len, 1);
        C(:,:,ff) = HmHm == HmHm';  % original def used in Brunet 2004
%         [rho,~] = corr(H);
%         C(:,:,ff) = corr(H);
        
%         figure;
%         imagesc(rho);
    end
    RE_max = max(RE_file);
    RE_min = min(RE_file);
    
    C_weighted = nan(size(C));
    for ff = 1:length(fname_all)
        C_weighted(:,:,ff) = C(:,:,ff)*(RE_max-RE_file(ff))/(RE_max-RE_min);
    end
    
    consensus = mean(C, 3);
    consensus_weighted = mean(C_weighted, 3);
    
    % Lower diagonal elements of consensus
    avec = [];
    avec_weighted = [];
    obs_len = size(consensus,1);
    for obs_i = 1:obs_len
        avec = [avec, consensus(obs_i, (obs_i+1):obs_len)];
        avec_weighted = [avec_weighted, consensus_weighted(obs_i, (obs_i+1):obs_len)];
    end
    
    % consensus entries are similarities, conversion to distances
    Y = 1-avec;
    Y_weighted = 1-avec_weighted;
    
    % Compute cophenetic coefficient
    Z = linkage(Y, 'average');
    [c,D] = cophenet(Z,Y);
    c_rank = [c_rank, c];
    Z_weighted = linkage(Y_weighted, 'average');
    [c_weighted,D_weighted] = cophenet(Z_weighted,Y_weighted);
    c_rank_weighted = [c_rank_weighted, c_weighted];
    
end

