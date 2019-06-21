% 2019 06 20

rank = 10;
beta = 0.1;
sm = 1e7;
lambda_all = [2,5,10,20];

for ilambda = 1:length(lambda_all)

    lambda = lambda_all(ilambda);
    A_str = sprintf('%d',lambda);

    filepath = '/media/wu-jung/internal_2tb/nmf_results/decomp_lambda_rank_sweep_revisit_1e4';    
    fname = sprintf('%s_r%02d_betaW%2.2f_betaH%2.2f_smoothness%0.2e_sparsity%0.2e.mat',...
                    'lambda_20190617_rank', rank, beta, beta, sm, lambda);
    A = load(fullfile(filepath,fname));

    step_str = {'1e2','2e2','5e2','1e3','2e3','5e3','1e4'};

    for istep = 1:length(step_str)
        step = eval(step_str{istep});

        % H at different steps
        H = squeeze(A.H_steps(:,:,step));
        for icomp=1:rank
            H_lines(icomp,:,istep) = H(icomp,:);
        end
        
        % W at different steps
        fig = figure;
        W = squeeze(A.W_steps(:,:,step));
        V = zeros(37,144*3,3);
        for icomp = 1:rank
            V = reshape(W(:,icomp),37,144*3);
            V = V/max(max(V));
            if rank>=8
                subplot(rank/2,2,icomp)
            else
                subplot(rank,1,icomp)
            end
            imagesc(V)
            caxis([0,1])
            colorbar
            title(sprintf('comp #%02d',icomp))
        end
        suptitle(['lambda = ',num2str(lambda),...
                  ', step = ',num2str(step)])
        saveas(gcf,fullfile(filepath,...
                            sprintf('lambda%s_W_step%05d.png',A_str,step)));
    end

    % Plot H
    for icomp=1:rank
        figure
        plot(squeeze(H_lines(icomp,:,:)));
        title(['lambda = ',num2str(lambda)]);
        legend(step_str)
        saveas(gcf,fullfile(filepath,sprintf('lambda%s_H_comp%d.png',A_str,icomp)));
    end

    % Plot objective
    figure
    plot(A.objective);
    ylabel('Objective');
    xlabel('Steps')
    title(['lambda = ',num2str(lambda)]);
    grid
    saveas(gcf,fullfile(filepath,sprintf('lambda%s_objective.fig',A_str)));
    saveas(gcf,fullfile(filepath,sprintf('lambda%s_objective.png',A_str)));

end

