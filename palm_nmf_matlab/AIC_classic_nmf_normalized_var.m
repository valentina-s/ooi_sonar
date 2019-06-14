% 2019 05 02  Test NMF with normalized variance and try to compute AIC



% Set paths
if ismac && isunix  % if on Mac
    addpath /Users/wu-jung/code_git/ooi_sonar/palm_nmf_matlab
    data_path = '/Users/wu-jung/code_git/ooi_sonar/sample_data/';
    save_path = '/Users/wu-jung/code_git/ooi_sonar/decomp_results_nmf_normvar/';
elseif isunix  % if on Linux
    addpath ~/code_git/ooi_sonar/palm_nmf_matlab
    data_path = '~/code_git/ooi_sonar/sample_data/';
    save_path = '~/code_git/ooi_sonar/decomp_results_nmf_normvar/';
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


% Run classic NMF, sweep through different rank
rep_all = 10;
rank_all = 1:30;
opt = statset('MaxIter',2000,'Display','final','TolFun',1e-6);

for rep = 1:rep_all
parfor rr = 1:length(rank_all)
    fprintf('rank = %d\n', rank_all(rr));
    
    % Original data
    [W, H, D] = nnmf(LL_norm,rank_all(rr),'option',opt,'algorithm','mult');
    save_file = sprintf('%s_r%02d_rep%04d.mat', ...
        'nmf_normalized_var_AIC', rank_all(rr), rep);
    save_file = fullfile(save_path, save_file);
    m = matfile(save_file,'writable',true);
    m.W = W;
    m.H = H;
    m.D = D;
end
end


% Compute AIC vs rank
AIC_rr = cell(length(rank_all),1);
for rr = 1:length(rank_all)
    rep_tmp = nan(length(rep_all),1);
    for rep = 1:rep_all
        save_file = sprintf('%s_r%02d_rep%04d.mat', ...
            'nmf_normalized_var_AIC', rank_all(rr), rep);

        % Load NMF results
        D = load(fullfile(save_path, save_file));

        % Calculate AIC
        fac1(rep) = norm(LL_norm - D.W*D.H, 'fro').^2;
        fac2(rep) = 2*rank_all(rr)*(size(LL,1)+size(LL,2));
        rep_tmp(rep) = fac1(rep) + fac2(rep);
    end
    AIC_rr{rr} = rep_tmp;
    AIC_rr_mean(rr) = mean(rep_tmp);
    AIC_rr_std(rr) = std(rep_tmp);
end


% Inspect components for a particular rank
rr = 3;
save_file = sprintf('%s_r%02d.mat', ...
    'nmf_normalized_var_AIC', rank_all(rr));
D = load(fullfile(save_path, save_file));

figure
for icomp=1:rr
    subplot(rr,1,icomp)
    imagesc(flipud(reshape(sigma_all.*D.W(:,icomp),37,144*3)));
end

figure
plot(D.H','linewidth',2)


