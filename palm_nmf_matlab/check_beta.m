% 2019 06 14

filepath = '/media/wu-jung/internal_2tb/nmf_results/decomp_beta_revisit_1e4/';
smoothness = 1e8;
A_str = '1e8';

switch smoothness
    case 1e2        
        A = load(fullfile(filepath,'beta_20190616_r03_betaW0.00_betaH0.00_smoothness000100.00.mat'));
    case 1e4
        A = load(fullfile(filepath,'beta_20190616_r03_betaW0.00_betaH0.00_smoothness010000.00.mat'));
    case 1e6
        A = load(fullfile(filepath,'beta_20190616_r03_betaW0.00_betaH0.00_smoothness1000000.00.mat'));
    case 1e8
        A = load(fullfile(filepath,'beta_20190616_r03_betaW0.00_betaH0.00_smoothness100000000.00.mat'));
end

comps = [1,2,3];
step_str = {'1e2','2e2','5e2','1e3','2e3','1e4'};

for istep = 1:length(step_str)
    step = eval(step_str{istep});
    step_idx = find(step == A.iter_save);
    
    % H at different steps
    H = squeeze(A.H_steps(:,:,step_idx));
    for icomp=1:length(comps)
        H_lines(icomp,:,istep) = H(icomp,:);
    end
    
    % W at different steps
    fig = figure;
    W = squeeze(A.W_steps(:,:,step_idx));
    V = zeros(37,144*3,3);
    for icomp = 1:length(comps)
        V = reshape(W(:,icomp),37,144*3);
        V = V/max(max(V));
        subplot(length(comps),1,icomp)
        imagesc(V)
        caxis([0,1])
        colorbar
        if icomp==1
            title(['smoothness = ',num2str(smoothness),...
                   ', step = ',num2str(step),...
                   ', betaH=betaW=0'])
        end
    end
    saveas(gcf,fullfile(filepath,...
                        sprintf('A%s_W_step%05d.png',A_str,step)));
end

% Plot H
for icomp=1:3
    figure
    plot(squeeze(H_lines(icomp,:,:)));
    title(sprintf('smoothness%s, betaH=betaW=0, comp #%d',...
                  A_str, icomp));
    legend(step_str)
    saveas(gcf,fullfile(filepath,sprintf('A%s_H_comp%d.png',A_str,icomp)));
end

% Plot objective
figure
plot(A.objective);
ylabel('Objective');
xlabel('Steps')
title(sprintf('smoothness%s, betaH=betaW=0',A_str));
if strcmp(A_str,'1e2')
    ylim([1.75,1.85]*1e6)
elseif strcmp(A_str,'1e4')
    ylim([1.5 2]*1e6)
elseif strcmp(A_str,'1e6')
    ylim([1,3]*1e6)
elseif strcmp(A_str,'1e8')
    ylim([4,8]*1e6)
end
grid
saveas(gcf,fullfile(filepath,sprintf('A%s_objective.fig',A_str)));
saveas(gcf,fullfile(filepath,sprintf('A%s_objective.png',A_str)));

