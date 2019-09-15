% Attempt to come up with a measure to select sparsity parameter


% Set paths
mvbs_path = '/Users/wu-jung/code_git/ooi_sonar/sample_data/';
mvbs_file = '20150817-20151017_MVBS_PCPcleaned.h5';
data_path1 = '/Volumes/MURI_4TB/nmf_results/ssNMF_sweep_sm_sp_20190623';
data_path2 = '/Volumes/MURI_4TB/nmf_results/ssNMF_sweep_sm_sp_20190626';
save_path = '/Volumes/MURI_4TB/nmf_results/ssNMF_sweep_sm_sp_check';


% Load MVBS data
curr_path = pwd;
cd(mvbs_path)  % has to switch folder to ensure reading h5 files correctly
L = h5read(mvbs_file,'/L');  % Load PCP-cleaned data
L = L-min(L(:));  % make it non-negative
cd(curr_path)


% Set params
r = 3;
betaH = 0.1;
betaW = 0.1;
max_iter = 2e4;
sp_all = [1,2,5,10,20,50];

% smoothness
% sm_order = [1,2,3,4];   % in folder ssNMF_sweep_sm_sp_20190626
% sm_order = [5,6,7,8];   % in folder ssNMF_sweep_sm_sp_20190623
sm_order = 1:8;
sm_all = repmat([1,2,5],length(sm_order),1);
for iorder = 1:length(sm_order)
    sm_all(iorder,:) = sm_all(iorder,:)*10^sm_order(iorder);
end
sm_all = sm_all';
sm_all = sm_all(:);


% Load final H
e_recon = nan(length(sm_all),length(sp_all),length(A.iter_save));
min_loc = nan(length(sm_all),length(sp_all));

for ism = 1:length(sm_all)

    sm = sm_all(ism);
    
    for isp = 1:length(sp_all)
        
        sp = sp_all(isp);
        
        fname = sprintf('%s_r%02d_betaW%2.2f_betaH%2.2f_smoothness%0.2e_sparsity%0.2e_maxiter%0.2e.mat', ...
            'ssNMF_runner', r, betaW, betaH, sm, sp, max_iter);
        if sm>5e4
            A = load(fullfile(data_path1,fname));  % sm_order = 5,6,7,8
        else
            A = load(fullfile(data_path2,fname));  % sm_order = 1,2,3,4
        end
        
        % get reconstruction error at max and last step
        for istep = 1:length(A.iter_save)
            e_recon(ism,isp,istep) = norm(L - A.W_steps(:,:,istep)*A.H_steps(:,:,istep), 'fro').^2;
        end
        
        % check location of objective minimum
        [mm, min_loc_tmp] = min(A.objective);
        min_loc(ism,isp) = min_loc_tmp;
    end
end

% Fig to check min_loc as a function of sp and sm
figure
imagesc(min_loc);
xticks(1:length(sp_all))
xticklabels(num2str(sp_all'))
yticks(1:length(sm_all))
yticklabels(num2str(sm_all))
xlabel('Sparsity param','fontsize',14)
ylabel('Smoothness param','fontsize',14)
tt = title('Step location of minimum objective');
set(tt,'fontsize',14)
colorbar
caxis([0,2e4])


