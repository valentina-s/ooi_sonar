% 2019 06 14

smoothness = 1e6;
A_str = '1e6';

switch smoothness
    case 1e2        
        A = load('smoothness_search_20190611_r03_betaW0.10_betaH0.10_smoothness000100.00.mat');
    case 1e4
        A = load('smoothness_search_20190611_r03_betaW0.10_betaH0.10_smoothness010000.00.mat');
    case 1e6
        A = load('smoothness_search_20190611_r03_betaW0.10_betaH0.10_smoothness1000000.00.mat');
    case 1e8
        A = load('smoothness_search_20190611_r03_betaW0.10_betaH0.10_smoothness100000000.00.mat');
end

comps = [1,2,3];
step_str = '2e3';
step = eval(step_str);

% W at different steps
fig = figure;
W = squeeze(A.W_steps(:,:,steps_all(istep)));
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
            ', step = ',num2str(steps_all(istep))])
    end
end
saveas(gcf,[A_str,'W_step',step_str,'.png']);

