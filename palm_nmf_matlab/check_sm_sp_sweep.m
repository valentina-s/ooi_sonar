% 2019 06 20

rank_all = 3:10;
betaH = 0.1;
betaW = 0.1;
max_iter = 2e4;
sp_all = [2,5,10,20];

% smoothness
sm_order = [6,7,8];
sm_all = repmat([1,2,5],length(sm_order),1);
for iorder = 1:length(sm_order)
    sm_all(iorder,:) = sm_all(iorder,:)*10^sm_order(iorder);
end
sm_all = sm_all';
sm_all = sm_all(:);
sm_all(2) = [];

step_str = {'1e2','2e2','5e2','1e3','2e3','5e3','1e4','2e4'};

data_path = '/Volumes/MURI_4TB/nmf_results/ssNMF_sweep_sm_sp_20190623';
save_path = '/Volumes/MURI_4TB/nmf_results/ssNMF_sweep_sm_sp_20190623/sm_sp_sweep';

% If save_path does not exist, create it
if ~exist(save_path, 'dir')
    mkdir(save_path)
end


for ism = 1:length(sm_all)
    for r = rank_all
        
        for isp = 1:length(sp_all)
            
            sm = sm_all(ism);
            sp = sp_all(isp);
            
            fname = sprintf('%s_r%02d_betaW%2.2f_betaH%2.2f_smoothness%0.2e_sparsity%0.2e_maxiter%0.2e.mat', ...
                'ssNMF_runner', r, betaW, betaH, sm, sp, max_iter);
            A = load(fullfile(data_path,fname));
            
            H_lines = zeros(r, size(A.H,2), length(step_str));
            for istep = 1:length(step_str)
                step = eval(step_str{istep});
                step_idx = find(step==A.iter_save);
                
                % H at different steps
                H = squeeze(A.H_steps(:,:,step_idx));
                for icomp=1:r
                    H_lines(icomp,:,istep) = H(icomp,:);
                end
                
                % W at different steps
                fig = figure;
                W = squeeze(A.W_steps(:,:,step_idx));
                V = zeros(37,144*3,3);
                for icomp = 1:r
                    V = reshape(W(:,icomp),37,144*3);
                    V = V/max(max(V));
                    if r>=8
                        subplot(ceil(r/2),2,icomp)
                    else
                        subplot(r,1,icomp)
                    end
                    imagesc(V)
                    caxis([0,1])
                    colorbar
                    title(sprintf('comp #%02d',icomp))
                end
                sgtitle(sprintf('r=%02d, sm=%.2e, sp=%.2e, step=%d',r,sm,sp,step));
                saveas(gcf,...
                    fullfile(save_path,...
                    sprintf('r%02d_sm%09.2e_sp%09.2e_W_step%05d.png',r,sm,sp,step)));
                close
            end
            
            % Plot H
            for icomp=1:r
                figure
                plot(squeeze(H_lines(icomp,:,:)));
                title(sprintf('comp#%02d, r=%02d, sm=%.2e, sp=%.2e',icomp,r,sm,sp));
                legend(step_str)
                saveas(gcf,...
                    fullfile(save_path,...
                    sprintf('r%02d_sm%09.2e_sp%09.2e_H_comp%02d.png',r,sm,sp,icomp)));
                close
            end
            
            % Plot objective
            figure
            plot(A.objective);
            ylabel('Objective');
            xlabel('Steps')
            title(sprintf('r=%02d, sm=%.2e, sp=%.2e',r,sm,sp));
            grid
            saveas(gcf,...
                fullfile(save_path,...
                sprintf('r%02d_sm%09.2e_sp%09.2e_objective.png',r,sm,sp)));
            saveas(gcf,...
                fullfile(save_path,...
                sprintf('r%02d_sm%09.2e_sp%09.2e_objective.fig',r,sm,sp)));
            close
        end
        
    end
end
