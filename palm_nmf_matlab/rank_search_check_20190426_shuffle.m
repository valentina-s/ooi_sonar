% 2019 04 26  Check param search results


if ismac && isunix  % if on Mac
    addpath /Users/wu-jung/code_git/ooi_sonar/palm_nmf_matlab/
    data_path = '/Users/wu-jung/code_git/ooi_sonar/sample_data/';
    save_path = '/Users/wu-jung/code_git/ooi_sonar/decomp_results/';
elseif isunix  % if on Linux
    addpath ~/code_git/ooi_sonar/palm_nmf_matlab/
    data_path = '~/code_git/ooi_sonar/sample_data/';
    save_path = '~/code_git/ooi_sonar/decomp_results/';
end

% Load PCP-cleaned data
data_file = '20150817-20151017_MVBS_PCPcleaned.h5';
cd(data_path)  % has to switch folder to ensure reading h5 files correctly
L = h5read(data_file,'/L');
L_sep = h5read(data_file,'/L_sep');
L_plot = h5read(data_file,'/L_plot');
LL = L-min(L(:));  % make it non-negative


% Shuffle observation
shuffle_opt = 1;

rng(168)
rand_seq = nan(size(LL));
LL_perm = nan(size(LL));
if shuffle_opt==0  % shuffle within row
    for ii=1:size(LL,1)
        rand_seq(ii,:) = randperm(size(LL,2));
        LL_perm(ii,:) = LL(ii,rand_seq(ii,:));
    end
elseif shuffle_opt==1  % shuffle within column
    for ii=1:size(LL,2)  
        rand_seq(:,ii) = randperm(size(LL,1));
        LL_perm(:,ii) = LL(rand_seq(:,ii),ii);
    end
end


% SVD
[U,S,V] = svd(LL);
[U_shuffle,S_shuffle,V_shuffle] = svd(LL_perm);

RMS_svd = nan(62,1);
RMS_svd_shuffle = nan(62,1);
for irank = 1:62
    RMS_svd(irank) = norm(LL - U(:,1:irank)*S(1:irank,1:irank)*V(:,1:irank)','fro');
    RMS_svd_shuffle(irank) = norm(LL_perm - U_shuffle(:,1:irank)*S_shuffle(1:irank,1:irank)*V_shuffle(:,1:irank)','fro');
end


% Sweep params
smoothness = 1000;
rank_all = 1:30;

objective_min = nan(length(rank_all),1);
objective_min_shuffle = nan(length(rank_all),1);
RMS_nmf = nan(length(rank_all),1);
RMS_nmf_shuffle = nan(length(rank_all),1);
for rr = 1:length(rank_all)
    % original data
    fname = sprintf('%s_r%02d_betaW%2.2f_betaH%2.2f_smoothness%08.2f.mat', ...
        'rank_search_20190426', rank_all(rr), ...
        0.1, 0.1, smoothness);
    D_original = load(fullfile(save_path, fname));
    objective_min(rr) = D_original.objective(end);
    RMS_nmf(rr) = norm(LL - D_original.W*D_original.H, 'fro');
    
    % shuffled data
    fname = sprintf('%s_r%02d_betaW%2.2f_betaH%2.2f_smoothness%08.2f.mat', ...
        'rank_search_20190426_shuffle', rank_all(rr), ...
        0.1, 0.1, smoothness);
    D_shuffle = load(fullfile(save_path, fname));
    objective_min_shuffle(rr) = D_shuffle.objective(end);
    LL_shuffle = nan(size(LL));
    for icol = 1:size(LL_shuffle,2)
        LL_shuffle(:,icol) = LL(D_shuffle.rand_seq(:,icol),icol);
    end
    RMS_nmf_shuffle(rr) = norm(LL_shuffle - D_shuffle.W*D_shuffle.H, 'fro');
end

% Plot
% RMS vs rank
figure
plot(1:62,RMS_svd,'o-')
hold on
plot(1:62,RMS_svd_shuffle,'o-')
plot(rank_all,RMS_nmf,'x-')
plot(rank_all,RMS_nmf_shuffle,'x-')
set(gca,'fontsize',14)
xlabel('Rank')
ylabel('RMS')
if shuffle_opt==0
    title('Shuffled within row','fontsize',16)
elseif shuffle_opt==1
    title('Shuffled within column','fontsize',16)
end
legend({'SVD','SVD shuffled','SNMF','SNMF shuffled'})
grid

% RMS slope vs rank
figure
plot(2:62, diff(RMS_svd), 'o-');
hold on
plot(2:62, diff(RMS_svd_shuffle), 'o-');
plot(rank_all(2:end),diff(RMS_nmf),'x-')
plot(rank_all(2:end),diff(RMS_nmf_shuffle),'x-')
set(gca,'fontsize',14)
xlabel('Rank')
ylabel('RMS')
if shuffle_opt==0
    title('Slope, shuffled within row','fontsize',16)
elseif shuffle_opt==1
    title('Slope, shuffled within column','fontsize',16)
end
legend({'SVD','SVD shuffled','SNMF','SNMF shuffled'},'location','best')
grid
