% 2019 06 14

filepath = '/media/wu-jung/internal_2tb/nmf_results/decomp_smoothness_revisit_2e4/';
smoothness = 1e8;
A_str = '1e8';

switch smoothness
    case 1e2        
        A = load(fullfile(filepath,'smoothness_search_20190611_r03_betaW0.10_betaH0.10_smoothness000100.00.mat'));
    case 1e4
        A = load(fullfile(filepath,'smoothness_search_20190611_r03_betaW0.10_betaH0.10_smoothness010000.00.mat'));
    case 1e6
        A = load(fullfile(filepath,'smoothness_search_20190611_r03_betaW0.10_betaH0.10_smoothness1000000.00.mat'));
    case 1e8
        A = load(fullfile(filepath,'smoothness_search_20190611_r03_betaW0.10_betaH0.10_smoothness100000000.00.mat'));
end

comps = [1,2,3];
step_str = {'1e2','2e2','1e3','2e3','1e4','2e4'};

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
            title(['smoothness = ',num2str(smoothness),...
                   ', step = ',num2str(step)])
        end
    end
    saveas(gcf,fullfile(filepath,...
                        sprintf('A%s_W_step%05d.png',A_str,step)));
end

% Plot H
for icomp=1:3
    figure
    plot(squeeze(H_lines(icomp,:,:)));
    legend(step_str)
    saveas(gcf,fullfile(filepath,sprintf('A%s_H_comp%d.png',A_str,icomp)));
end

% Plot objective
figure
plot(A.objective);
ylabel('Objective');
xlabel('Steps')
if strcmp(A_str,'1e2')
    ylim([1.75,1.85]*1e6)
elseif strcmp(A_str,'1e4')
    ylim([1.8 2]*1e6)
elseif strcmp(A_str,'1e6')
    ylim([0.23,0.3]*1e7)
elseif strcmp(A_str,'1e8')
    ylim([5,6]*1e6)
end
grid
saveas(gcf,fullfile(filepath,sprintf('A%s_objective.png',A_str)));

