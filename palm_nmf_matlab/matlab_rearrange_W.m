% 2017 11 28  Rename files to form better sequence

clear

if ismac
    data_path = '/Users/wujung/Downloads/';
    save_path = '/Users/wujung/Downloads/';
    addpath /Users/wujung/code/ooi_sonar/
else
    data_path = '/media/wu-jung/wjlee_apl_2/ooi_zplsc_palm_nmf_2018/';
    save_path = '/media/wu-jung/wjlee_apl_2/ooi_zplsc_palm_nmf_2018/';
    addpath ~/internal_2tb/ooi_sonar/ooi_sonar_code/
end

files = dir(fullfile(data_path,'palm_nmf_r3_*.mat'));

A = load(fullfile(data_path,files(1).name));

n_comp = size(A.W,2);

idx_match = nan(n_comp,2);
idx_match_all = nan(length(files),n_comp);
fig = figure('position',[100,200,500,600]);
for iF=1:length(files)
    A = load(fullfile(data_path,files(iF).name));

    s = strsplit(files(iF).name,'.mat');
    
    % distinguish if beta or eta param first in filename
    if strcmp(s{1}(13),'e')
        eta = str2double(s{1}(16:end));
    else
        betaW = str2double(s{1}(18:20));
        betaH = str2double(s{1}(27:29));
        eta = str2double(s{1}(34:end));
    end

    if iF~=1
        RHO = corr(W_prev,A.W);
        for iC=1:n_comp
            [~,idx] = max(RHO(:));
            [idx_match(iC,1),idx_match(iC,2)] = ind2sub(size(RHO),idx);
            RHO(idx_match(iC,1),:) = -Inf;
            RHO(:,idx_match(iC,2)) = -Inf;
        end
        idx_match = sortrows(idx_match);
    else
        idx_match(:,1) = [1:n_comp]';
        idx_match(:,2) = [1:n_comp]';
    end
    idx_match_all(iF,:) = idx_match(:,2);
    W_prev = A.W(:,idx_match(:,2));

    W1_all(:,iF) = W_prev(:,1);
    W2_all(:,iF) = W_prev(:,2);
    W3_all(:,iF) = W_prev(:,3);
    if n_comp==4
        W4_all(:,iF) = W_prev(:,4);
    end

    V = reshape(W_prev,37,144*3,n_comp);

    figure(fig);
%     figure;
    for iC=1:n_comp
        subplot(n_comp+2,1,iC)
        imagesc(V(:,:,iC))
        axis xy
        colorbar
        if iC==1
            title(['eta=',num2str(eta)])
        end
        hold on
        plot([144,144],[0,37],'w--','linewidth',1);
        plot([288,288],[0,37],'w--','linewidth',1);
    end
    subplot(n_comp+2,1,n_comp+(1:2))
    plot(A.H(idx_match(:,2),:)','linewidth',2);
    legend(num2str([1:n_comp]'))

    saveas(gcf,fullfile(save_path,sprintf('%s.png',s{1})),'png');
end


