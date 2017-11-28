% 2017 11 27  Test PALM-NMF algorithm on MVBS data

clear

data_path = '/Users/wujung/Downloads/';
data_file = '20150817-20151017_MVBS_PCPcleaned.h5';

cd(data_path)  % has to switch folder to ensure reading h5 files correctly
%h5disp(data_file)


L = h5read(data_file,'/L');
L_sep = h5read(data_file,'/L_sep');
L_plot = h5read(data_file,'/L_plot');
% S = h5read(data_file,'/S');
% S_sep = h5read(data_file,'/S_sep');
% S_plot = h5read(data_file,'/S_plot');
depth_bin_size = h5read(data_file,'/depth_bin_size');
ping_time = h5read(data_file,'/ping_time');
ping_per_day_mvbs = h5read(data_file,'/ping_per_day_mvbs');

params.r = 4;
params.max_iter = 100;
[W_100, H_100, objective_100, iter_times_100] = palm_nmf(L-min(L(:)), params);

params.r = 4;
params.max_iter = 100;
params.betaH = 1;
params.betaW = 0.1;
[W_100_betaH1, H_100_betaH1, objective_100_betaH1, iter_times_100_betaH1] = palm_nmf(L-min(L(:)), params);

% params.r = 4;
% params.max_iter = 500;
% [W_500, H_500, objective_500, iter_times_500] = palm_nmf(L-min(L(:)), params);

params.r = 4;
params.max_iter = 100;
params.betaH = 1;
params.betaW = 1;
[W_100_betaH1_betaW1, H_100_betaH1_betaW1, objective_100_betaH1_betaW1, iter_times_100_betaH1_betaW1] = palm_nmf(L-min(L(:)), params);


V_100 = reshape(W_100,37,144*3,[],4);
% V_500 = reshape(W_500,37,144*3,[],4);
V_100_betaH1 = reshape(W_100_betaH1,37,144*3,[],4);
V_100_betaH1_betaW1 = reshape(W_100_betaH1_betaW1,37,144*3,[],4);


figure;
for iC=1:4
    subplot(4,1,iC)
    imagesc(V_100(:,:,iC))
    axis xy
    colorbar
    if iC==1
        title('100 iter, betaH=betaW=0.1')
    end
    hold on
    plot([144,144],[0,37],'w--','linewidth',1);
    plot([288,288],[0,37],'w--','linewidth',1);
end


figure;
for iC=1:4
    subplot(4,1,iC)
    imagesc(V_100_betaH1(:,:,iC))
    axis xy
    colorbar
    if iC==1
        title('100 iter, betaH=1, betaW=0.1')
    end
    hold on
    plot([144,144],[0,37],'w--','linewidth',1);
    plot([288,288],[0,37],'w--','linewidth',1);
end


figure;
for iC=1:4
    subplot(4,1,iC)
    imagesc(V_100_betaH1_betaW1(:,:,iC))
    axis xy
    colorbar
    if iC==1
        title('100 iter, betaH=1, betaW=1')
    end
    hold on
    plot([144,144],[0,37],'w--','linewidth',1);
    plot([288,288],[0,37],'w--','linewidth',1);
end





