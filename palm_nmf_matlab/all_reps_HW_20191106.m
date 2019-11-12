% Plot all H and W results from the same param combination

curr_path = pwd;


% Get normalized vector
sv_path = '../sample_data';
sv_file = '20150817-20151017_MVBS_PCPcleaned.h5';
cd(sv_path)  % has to switch folder to ensure reading h5 files correctly
L = h5read(sv_file,'/L');
sigma_all = std(L-min(L(:)),0,2);  % variance of non-negative values
cd(curr_path)


% Paths
data_path = '/Volumes/MURI_4TB/nmf_results/ssNMF_20191106/';
save_path = '/Volumes/MURI_4TB/nmf_results/all_reps_HW_20191106';


% ss-NMF parameters
betaHW = 0.1;
sm = 5e4;
max_iter = 2e4;


% If save_path does not exist, create it
if ~exist(save_path, 'dir')
    mkdir(save_path)
end


% Gather H and W with the smallest objective
rep_num = 10;
pc_per_fig = 3;  % param combinations to plot per figure
pc_cnt = 1;      % flag to control if to initialize new figure

for r = 4
for sp = [0.1, 0.2, 0.5, 1, 2, 5]
    
sg_str = sprintf('rank=%2d, sm=%.2e, sp=%.1f',...
    r, sm, sp);    

% If save_path does not exist, create it
if ~exist(save_path, 'dir')
    mkdir(save_path)
end

fig_str = sprintf('r%02d_betaHW%2.2d_sm%0.2e_sp%0.2e', r, betaHW, sm, sp);


tic

for irep = 0+(1:rep_num)
    
    % Assemble file name
    fname = sprintf('%s_r%02d_betaW%2.2f_betaH%2.2f_sm%0.2e_sp%0.2e_maxiter%0.2e_rep%03d.mat', ...
        'ssNMF_runner_normvar_repeat', r, betaHW, betaHW, sm, sp, max_iter, irep);

    % Load data
    disp(['Loading ', fname]);
    A = matfile(fullfile(data_path,fname));  % enable dynamic loading
    
    [min_obj_at_step, min_obj_step] = min(A.objective);
    [~, min_obj_idx] = min(min_obj_step-A.iter_save);  % index of the closest step saved
    H_min = squeeze(A.H_steps(:,:,min_obj_idx));
    W_min = squeeze(A.W_steps(:,:,min_obj_idx));
    
    % Save to overall results
    S.min_obj_step(irep) = min_obj_step;                         % the actual step with smallest objective
    S.min_obj_step_closest(irep) = A.iter_save(min_obj_idx,1);   % the closest step saved
    S.min_obj_at_step(irep) = min_obj_at_step;                   % the actual smallest objective
    S.min_obj_at_step_closest(irep) = ...
        A.objective(A.iter_save(min_obj_idx,1),1);               % the saved smallest objective
    S.H_min(:,:,irep) = H_min;
    S.W_min(:,:,irep) = W_min;
end
S.sigma_all = sigma_all;   % variance of each pixel

toc


% Find matching between components across reps
pair_all = zeros(r,rep_num);
pair_all(:,1) = 1:r;
for irep = 2:rep_num
    H_min_rep = squeeze(S.H_min(:,:,1));
    H_min_rep_norm = H_min_rep/max(H_min_rep);
    hcorr = corr(squeeze(S.H_min(:,:,1))', ...  % correlation between columns
                 squeeze(S.H_min(:,:,irep))');

    pair = zeros(r,2);  % matrix to store corresponding pair of components
    for ir = 1:r  % go through each component
        [~, mm_idx] = max(hcorr(:));
        [pair(ir,1),pair(ir,2)] = ind2sub([r, r], mm_idx);
        hcorr(pair(ir,1),:) = 0;
    end
    pair = sortrows(pair,1);  % sort rows based on first column
    
    pair_all(:,irep) = (pair(:,2));
end


% Check re-ordered H
fig_h_all = figure;
set(0,'DefaultAxesColorOrder',brewermap(rep_num,'Blues'))
hold on
ih = 4;
for irep = 1:rep_num
    plot(S.H_min(pair_all(ih,irep),:,irep),'linewidth',1.5);
end
set(0,'DefaultAxesColorOrder','remove')

% Find the rep that has the minimum objective
[~, rep_min_obj_idx] = min(S.min_obj_at_step_closest);
plot(S.H_min(pair_all(ih,rep_min_obj_idx),:,rep_min_obj_idx),'r','linewidth',2,'linestyle','--');


% ===================== WORK UP TO HERE =======================

% Plot W_min - normalized
fig_W_norm = figure;
V = zeros(37,144*3,3);
for icomp = 1:r
    V = reshape(squeeze(S.W_min(:,icomp,rep_min_obj_idx)),37,144*3);
    if r>=6  % make fig 2 columns for better viz
        subplot(ceil(r/2),2,icomp)
    else
        subplot(r,1,icomp)
    end
    imagesc(V)
    colorbar
    title(sprintf('comp #%02d',icomp))
end
xlabel('Hour of day x 3 freq','fontsize',12);
ylabel('Depth bins','fontsize',12);
sgtitle(['normalized, ', sg_str])


% Plot W_min - rescaled
fig_W_res = figure;
V = zeros(37,144*3,3);
for icomp = 1:r
    V = reshape(squeeze(S.W_min(:,icomp,rep_min_obj_idx)).*sigma_all,37,144*3);
    if r>=6  % make fig 2 columns for better viz
        subplot(ceil(r/2),2,icomp)
    else
        subplot(r,1,icomp)
    end
    imagesc(V)
    colorbar
    title(sprintf('comp #%02d',icomp))
end
xlabel('Hour of day x 3 freq','fontsize',12);
ylabel('Depth bins','fontsize',12);
sgtitle(['rescaled, ', sg_str])


% Plot H_min
fig_H = figure;
plot(squeeze(S.H_min(:,:,rep_min_obj_idx))', 'linewidth',2);
ll = legend();
set(ll,'fontsize',14)
xlabel('Day of observation','fontsize',12);
ylabel('Activation','fontsize',12);
sgtitle(sg_str)


% Plot objective distribution
fig_obj = figure;
subplot(211)
hist(S.min_obj_at_step_closest,30);
xlabel('Objective','fontsize',12);
ylabel('Count','fontsize',12);
subplot(212)
plot(S.min_obj_at_step_closest,'o-');
xlabel('Run number','fontsize',12);
ylabel('Min objective','fontsize',12);
title(sg_str)


% Save summary figs and and results
saveas(fig_W_norm, fullfile(save_path, [fig_str, '_Wmin_normalized.png']));
saveas(fig_W_norm, fullfile(save_path, [fig_str, '_Wmin_normalized.fig']));
saveas(fig_W_res, fullfile(save_path, [fig_str, '_Wmin_rescaled.png']));
saveas(fig_W_res, fullfile(save_path, [fig_str, '_Wmin_rescaled.fig']));
saveas(fig_H, fullfile(save_path, [fig_str, '_Hmin.png']));
saveas(fig_H, fullfile(save_path, [fig_str, '_Hmin.fig']));
saveas(fig_obj, fullfile(save_path, [fig_str, '_obj.png']));
saveas(fig_obj, fullfile(save_path, [fig_str, '_obj.fig']));
save(fullfile(save_path, [fig_str, '_summary.mat']), '-struct', 'S');


sprintf('Smallest objective: %.3f', rep_min_obj);

end
end
