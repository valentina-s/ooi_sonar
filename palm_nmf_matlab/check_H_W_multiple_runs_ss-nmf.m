% 2019 05 02  Check correlation across multiple NMF runs in terms of 
%             components and coefficients


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


% Check correlation between components
rank = 4;
rep_all = 1:10;

% Get results from first NMF
save_file = sprintf('rank_fixed_20190502_repeat_r%02d_betaW0.10_betaH0.10_smoothness100.00_rep0001.mat',rank);
A = load(fullfile(save_path, save_file));
W0 = A.W;
H0 = A.H;
D0 = A.objective(end);

seq = nan(length(rep_all),rank);
residual_all = nan(length(rep_all),1);
for rep = rep_all
    save_file = sprintf('%s_r%02d_betaW0.10_betaH0.10_smoothness100.00_rep%04d.mat', ...
        'rank_fixed_20190502_repeat', rank, rep);

    % Load NMF results
    A = load(fullfile(save_path, save_file));

    % Get re-order sequence
%     [rho,~] = corr(W0,A.W);  % re-order based on W
    [rho,~] = corr(H0',A.H');  % re-order based on H
    seq(rep,:) = find_match_factor_seq(rho,rank);
    
    % Error
    residual_all(rep) = A.objective(end);
end


% Plot multiple Ws for insepction
fig_W = figure;
fig_H = figure;
for rep = rep_all
    save_file = sprintf('%s_r%02d_betaW0.10_betaH0.10_smoothness100.00_rep%04d.mat', ...
        'rank_fixed_20190502_repeat', rank, rep);

    % Load NMF results
    A = load(fullfile(save_path, save_file));
    
    % Re-order W and H
    W_new = A.W(:,seq(rep,:));
    H_new = A.H(seq(rep,:),:);
    
    % Plot
    W_plot = [];
    for icomp = 1:rank
        W_plot = [W_plot;flipud(reshape(sigma_all.*W_new(:,icomp),37,144*3))];
    end
    figure(fig_W)
    subplot(rep_all(end)/2,2,rep)
    imagesc(W_plot);
    title(sprintf('rep = %02d',rep))
    
    figure(fig_H)
    subplot(rep_all(end)/2,2,rep)
    plot(H_new','linewidth',2)
    title(sprintf('rep = %02d',rep))
end






