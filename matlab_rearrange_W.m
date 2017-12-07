% 2017 11 28  Rename files to form better sequence

clear

if ismac
    data_path = '/Users/wujung/Downloads/';
    addpath /Users/wujung/code/ooi_sonar/
else
    data_path = '/media/wu-jung/wjlee_apl_2/ooi_zplsc_new/';
    addpath ~/internal_2tb/ooi_sonar/ooi_sonar_code/
end

files = dir(fullfile(data_path,'palm_nmf_r3_*.mat'));

A = load(fullfile(data_path,files(1).name));

idx_match = nan(size(A.W,2),2);
idx_match_all = nan(length(files),size(A.W,2));
fig = figure('position',[100,200,500,600]);
for iF=1:length(files)
    A = load(fullfile(data_path,files(iF).name));

    s = strsplit(files(iF).name,'.mat');
    eta = str2double(s{1}(16:end));

    if iF~=1
        RHO = corr(W_prev,A.W);
        for iC=1:size(A.W,2)
            [~,idx] = max(RHO(:));
            [idx_match(iC,1),idx_match(iC,2)] = ind2sub(size(RHO),idx);
            RHO(idx_match(iC,1),:) = -Inf;
            RHO(:,idx_match(iC,2)) = -Inf;
        end
        idx_match = sortrows(idx_match);
    else
        idx_match(:,1) = [1:size(A.W,2)]';
        idx_match(:,2) = [1:size(A.W,2)]';
    end
    idx_match_all(iF,:) = idx_match(:,2);
    W_prev = A.W(:,idx_match(:,2));

    W1_all(:,iF) = W_prev(:,1);
    W2_all(:,iF) = W_prev(:,2);
    W3_all(:,iF) = W_prev(:,3);
    if size(A.W,2)==4
        W4_all(:,iF) = W_prev(:,4);
    end

    V = reshape(W_prev,37,144*3,size(A.W,2));

    figure(fig);
%     figure;
    for iC=1:size(A.W,2)
        subplot(size(A.W,2)+2,1,iC)
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
    subplot(size(A.W,2)+2,1,size(A.W,2)+(1:2))
    plot(A.H(idx_match(:,2),:)','linewidth',2);
    legend('1','2','3','4');

    saveas(gcf,fullfile(data_path,sprintf('palm_nmf_eta%06.2f.png',eta)),'png');
end
