% 2019 06 14

sm = 1e6;
lambda_all = [1,2,5,10,20,50];

for ilambda = 1:length(lambda_all)

    lambda = lambda_all(ilambda);
    A_str = sprintf('%d',lambda);

    filepath = '/media/wu-jung/internal_2tb/nmf_results/decomp_lambda/';
    fname = sprintf('lambda_20190616_r03_betaW0.10_betaH0.10_smoothness%09.2f_sparsity%09.2f.mat',sm,lambda);
    A = load(fullfile(filepath,fname));

    comps = [1,2,3];
    step_str = {'1e2','2e2','5e2','1e3'};

    for istep = 1:length(step_str)
        step = eval(step_str{istep});

        % H at different steps
        H = squeeze(A.H_steps(:,:,step));
        for icomp=1:length(comps)
            H_lines(icomp,:,istep) = H(icomp,:);
        end
        
        % W at different steps
        fig = figure;
        W = squeeze(A.W_steps(:,:,step));
        V = zeros(37,144*3,3);
        for icomp = 1:length(comps)
            V = reshape(W(:,icomp),37,144*3);
            V = V/max(max(V));
            subplot(length(comps),1,icomp)
            imagesc(V)
            caxis([0,1])
            colorbar
            if icomp==1
                title(['lambda = ',num2str(lambda),...
                       ', step = ',num2str(step)])
            end
        end
        saveas(gcf,fullfile(filepath,...
                            sprintf('lambda%s_W_step%05d.png',A_str,step)));
    end

    % Plot H
    for icomp=1:3
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

