% Attempt to come up with a measure to select smoothness parameter


% Set paths
mvbs_path = '/Users/wu-jung/code_git/ooi_sonar/sample_data/';
mvbs_file = '20150817-20151017_MVBS_PCPcleaned.h5';
data_path = '/Volumes/MURI_4TB/nmf_results/ssNMF_sweep_sm_sp_20190623';
% data_path = '/Volumes/MURI_4TB/nmf_results/ssNMF_sweep_sm_sp_20190626';
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
sp = 20;

% smoothness
% sm_order = [1,2,3,4];   % in folder ssNMF_sweep_sm_sp_20190626
sm_order = [5,6,7,8];   % in folder ssNMF_sweep_sm_sp_20190623
sm_all = repmat([1,2,5],length(sm_order),1);
for iorder = 1:length(sm_order)
    sm_all(iorder,:) = sm_all(iorder,:)*10^sm_order(iorder);
end
sm_all = sm_all';
sm_all = sm_all(:);
sm_all(2) = [];


% Load final H
sm_max = nan(length(sm_all),1);
sm_last = nan(length(sm_all),1);
e_recon_max = nan(length(sm_all),1);
e_recon_last = nan(length(sm_all),1);

for ism = 1:length(sm_all)
    
    sm = sm_all(ism);
    
    fname = sprintf('%s_r%02d_betaW%2.2f_betaH%2.2f_smoothness%0.2e_sparsity%0.2e_maxiter%0.2e.mat', ...
        'ssNMF_runner', r, betaW, betaH, sm, sp, max_iter);
    A = load(fullfile(data_path,fname));
    
    % normalize within each H step
    A.H_steps_norm = nan(size(A.H_steps));
    for istep = 1:length(A.iter_save)
        for icomp = 1:A.params.r
            A.H_steps_norm(icomp,:,istep) = A.H_steps(icomp,:,istep)/max(A.H_steps(icomp,:,istep));
        end
    end
    
    % smoothenss measure
    A.H_sm_corr = nan([A.params.r,length(A.iter_save)]);
    for icomp = 1:A.params.r
        for istep = 1:length(A.iter_save)
            corr_tmp = xcorr(squeeze(A.H_steps_norm(icomp,:,istep)),1,'coeff');
            A.H_sm_corr(icomp,istep) = corr_tmp(1);
        end
    end
    
%     % check smoothness across steps and components
%     figure
%     plot(A.iter_save,A.H_sm_corr,'-x')
%     set(gca,'xscale','log')
%     xlabel('Iteration','fontsize',14)
%     ylabel('Smoothness','fontsize',14)

    % get averaged smoothness at max and last step
    mm = mean(A.H_sm_corr);
    [sm_max(ism), max_loc] = max(mm);
    sm_last(ism) = mm(end);
    
    % get reconstruction error at max and last step
    e_recon_max(ism) = norm(L - A.W_steps(:,:,max_loc)*A.H_steps(:,:,max_loc), 'fro').^2;
    e_recon_last(ism) = norm(L - A.W_steps(:,:,end)*A.H_steps(:,:,end), 'fro').^2;
end

figure('position',[672    53   560   652])
subplot(311)  % smoothness variation as a function of smoothness param
plot(sm_all,sm_max,'-x')
hold on
plot(sm_all,sm_last,'-x')
ylabel('Smoothness','fontsize',14)
grid
ll = legend('sm max xcorr lag-1','sm last xcorr lag-1');
set(ll,'fontsize',14,'location','best');
title('Smoothness variation')

subplot(312)  % 1-step difference on smoothness varitation curve
plot(sm_all(2:end),diff(sm_max),'-x')
hold on
plot(sm_all(2:end),diff(sm_last),'-x')
grid
ylabel('Difference','fontsize',14)
title('Diff of smoothness variation')

subplot(313)  % reconstruction error
plot(sm_all,e_recon_max,'-x')
hold on
plot(sm_all,e_recon_last,'-x')
grid
ylabel('Frobenius norm of diff','fontsize',14)
title('Reconstruction error')
xlabel('Smoothness parameter','fontsize',14)

tt = sgtitle(sprintf('Sparsity = %d', sp));
set(tt,'fontsize',14);


