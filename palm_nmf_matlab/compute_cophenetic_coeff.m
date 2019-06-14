% 2019 04 29  hierarchical clustering and cophenetic coeff for NMF results
%             over multiple runs


if ismac && isunix  % if on Mac
    addpath /Users/wu-jung/code_git/ooi_sonar/palm_nmf_matlab/
    data_path = '/Users/wu-jung/code_git/ooi_sonar/sample_data/';
    save_path = '/Users/wu-jung/code_git/ooi_sonar/decomp_results_repeat/';
elseif isunix  % if on Linux
    addpath ~/code_git/ooi_sonar/palm_nmf_matlab/
    data_path = '~/code_git/ooi_sonar/sample_data/';
    save_path = '~/code_git/ooi_sonar/decomp_results_repeat/';
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


% Load 1 file for params
fname = sprintf('%s_r%02d_betaW%2.2f_betaH%2.2f_smoothness%06.2f_rep%04d.mat', ...
    'rank_fixed_20190427_repeat', 3, 0.1, 0.1, 100, 1);
D = load(fullfile(save_path, fname));
obs_len = size(D.H,2);
clear D


% Params
smoothness = 100;
rank_all = 2:8;
c_rank = [];
c_rank_weighted = [];
for rank = rank_all
    
    % Get filenames of results with the same param combination
    fname_pre = sprintf('rank_fixed_20190427_repeat_r%02d_betaW0.10_betaH0.10_smoothness%06.2f_rep*.mat', ...
        rank, smoothness);
    fname_all = dir(fullfile(save_path,fname_pre));
    
    % Concensus matrix
    C = nan(obs_len, obs_len, length(fname_all));
    RE_file = nan(length(fname_all),1);
    for ff = 1:length(fname_all)
        % Load file
        fname = fname_all(ff).name;
        load(fullfile(save_path, fname));
        
        % Get RE
        RE_file(ff) = norm(LL - W*H, 'fro');
        
        % Connectivity matrix
        [~, Hm_idx] = max(H, [], 1);  % membership assignment
        HmHm = repmat(Hm_idx, obs_len, 1);
        C(:,:,ff) = HmHm == HmHm';
    end
    RE_max = max(RE_file);
    RE_min = min(RE_file);
    
    C_weighted = nan(size(C));
    for ff = 1:length(fname_all)
        C_weighted(:,:,ff) = C(:,:,ff)*(RE_max-RE_file(ff))/(RE_max-RE_min);
    end
    
    consensus = mean(C, 3);
    consensus_weighted = mean(C_weighted, 3);
    
%     figure
%     imagesc(consensus);
%     title(sprintf('Rank = %2d', rank), 'fontsize', 14);
%     colorbar
    
    % Lower diagonal elements of consensus
    avec = [];
    avec_weighted = [];
    obs_len = size(consensus,1);
    for obs_i = 1:obs_len
        avec = [avec, consensus((obs_i+1):obs_len, obs_i)'];
        avec_weighted = [avec_weighted, consensus_weighted((obs_i+1):obs_len, obs_i)'];
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