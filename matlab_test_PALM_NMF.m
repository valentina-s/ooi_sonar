% 2017 11 27  Test PALM-NMF algorithm on MVBS data

clear

addpath ~/internal_2tb/ooi_sonar/ooi_sonar_code/

data_path = '/media/wu-jung/wjlee_apl_2/ooi_zplsc_new/';
data_file = '20150817-20151017_MVBS_PCPcleaned.h5';

cd(data_path)  % has to switch folder to ensure reading h5 files correctly
%h5disp(data_file)

% Load PCP-cleaned data
L = h5read(data_file,'/L');
L_sep = h5read(data_file,'/L_sep');
L_plot = h5read(data_file,'/L_plot');
% S = h5read(data_file,'/S');
% S_sep = h5read(data_file,'/S_sep');
% S_plot = h5read(data_file,'/S_plot');
depth_bin_size = h5read(data_file,'/depth_bin_size');
ping_time = h5read(data_file,'/ping_time');
ping_per_day_mvbs = h5read(data_file,'/ping_per_day_mvbs');
depth_bin_num = 37;


% iteration over eta (temporal continuity while holding other coefs=0)
eta_all = [0:0.1:5,6:10,10:5:100]';
len = length(eta_all);
f1 = 'r';
v1 = mat2cell(4*ones(len,1),ones(len,1));
f2 = 'max_iter';
v2 = mat2cell(300*ones(len,1),ones(len,1));
f3 = 'betaW';
v3 = mat2cell(0*ones(len,1),ones(len,1));
f4 = 'betaH';
v4 = mat2cell(0*ones(len,1),ones(len,1));
f5 = 'eta';
v5 = mat2cell(eta_all,ones(len,1));
params_all = struct(f1,v1,f2,v2,f3,v3,f4,v4,f5,v5);

LL = L-min(L(:));

for iM = 1:25
    fprintf('eta=%2.1f\n',eta_all(iM))
    [W, H, objective, iter_times] = ...
        palm_nmf(LL, params_all(iM));
    W_all{iM} = W;
    H_all{iM} = H;
    objective_all{iM} = objective;
    params = params_all(iM);
    save(sprintf('palm_nmf_eta%2.1f.mat',eta_all(iM)),...
         'W','H','objective','iter_times','params');
end


if 0
figure;
for iC=1:4
    subplot(4,1,iC)
    imagesc(V_100_H1_W1(:,:,iC))
    axis xy
    colorbar
    if iC==1
        title('100 iter, betaH=1, betaW=1')
    end
    hold on
    plot([144,144],[0,37],'w--','linewidth',1);
    plot([288,288],[0,37],'w--','linewidth',1);
end
end




