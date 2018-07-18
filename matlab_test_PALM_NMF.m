% 2017 11 27  Test PALM-NMF algorithm on MVBS data

clear

addpath ~/internal_2tb/ooi_sonar/ooi_sonar_code/

data_path = '/media/wu-jung/wjlee_apl_2/ooi_zplsc_current/';
save_path = '/media/wu-jung/wjlee_apl_2/ooi_zplsc_palm_nmf/';
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


LL = L-min(L(:));  % make it non-negative


% Run PALM-NMF

% r = 3
eta_all = [0.01:0.01:0.09,0.1:0.1:1,10,100]';
len = length(eta_all);
f1 = 'r';
v1 = mat2cell(3*ones(len,1),ones(len,1));
f2 = 'max_iter';
v2 = mat2cell(300*ones(len,1),ones(len,1));
f3 = 'betaW';
v3 = mat2cell(0*ones(len,1),ones(len,1));
f4 = 'betaH';
v4 = mat2cell(0*ones(len,1),ones(len,1));
f5 = 'eta';
v5 = mat2cell(eta_all,ones(len,1));
params_all = struct(f1,v1,f2,v2,f3,v3,f4,v4,f5,v5);

for iM = 1:len
    fprintf('eta=%6.2f\n',eta_all(iM))
    [W, H, objective, iter_times] = ...
        palm_nmf(LL, params_all(iM));
    W_all{iM} = W;
    H_all{iM} = H;
    objective_all{iM} = objective;
    params = params_all(iM);
    save_file = sprintf('palm_nmf_r3_betaW0.0_betaH0.0_eta%06.2f.mat',eta_all(iM));
    save(fullfile(save_path,save_file),...
         'W','H','objective','iter_times','params');
end


eta_all = [0.01:0.01:0.09,0.1:0.1:1,10,100]';
len = length(eta_all);
f1 = 'r';
v1 = mat2cell(3*ones(len,1),ones(len,1));
f2 = 'max_iter';
v2 = mat2cell(300*ones(len,1),ones(len,1));
f3 = 'betaW';
v3 = mat2cell(0.1*ones(len,1),ones(len,1));
f4 = 'betaH';
v4 = mat2cell(0.1*ones(len,1),ones(len,1));
f5 = 'eta';
v5 = mat2cell(eta_all,ones(len,1));
params_all = struct(f1,v1,f2,v2,f3,v3,f4,v4,f5,v5);

for iM = 1:len
    fprintf('eta=%6.2f\n',eta_all(iM))
    [W, H, objective, iter_times] = ...
        palm_nmf(LL, params_all(iM));
    W_all{iM} = W;
    H_all{iM} = H;
    objective_all{iM} = objective;
    params = params_all(iM);
    save_file = sprintf('palm_nmf_r3_betaW0.1_betaH0.1_eta%06.2f.mat',eta_all(iM));
    save(fullfile(save_path,save_file),...
         'W','H','objective','iter_times','params');
end


if 0
% r = 4
eta_all = [0.01:0.01:0.09,0.1:0.1:1,10,100]';
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

for iM = 1:len
    fprintf('eta=%6.2f\n',eta_all(iM))
    [W, H, objective, iter_times] = ...
        palm_nmf(LL, params_all(iM));
    W_all{iM} = W;
    H_all{iM} = H;
    objective_all{iM} = objective;
    params = params_all(iM);
    save_file = sprintf('palm_nmf_r4_betaW0.0_betaH0.0_eta%06.2f.mat',eta_all(iM));
    save(fullfile(save_path,save_file),...
         'W','H','objective','iter_times','params');
end


eta_all = [0.01:0.01:0.09,0.1:0.1:1,10,100]';
len = length(eta_all);
f1 = 'r';
v1 = mat2cell(4*ones(len,1),ones(len,1));
f2 = 'max_iter';
v2 = mat2cell(300*ones(len,1),ones(len,1));
f3 = 'betaW';
v3 = mat2cell(0.1*ones(len,1),ones(len,1));
f4 = 'betaH';
v4 = mat2cell(0.1*ones(len,1),ones(len,1));
f5 = 'eta';
v5 = mat2cell(eta_all,ones(len,1));
params_all = struct(f1,v1,f2,v2,f3,v3,f4,v4,f5,v5);

for iM = 1:len
    fprintf('eta=%06.2f\n',eta_all(iM))
    [W, H, objective, iter_times] = ...
        palm_nmf(LL, params_all(iM));
    W_all{iM} = W;
    H_all{iM} = H;
    objective_all{iM} = objective;
    params = params_all(iM);
    save_file = sprintf('palm_nmf_r4_betaW0.1_betaH0.1_eta%06.2f.mat',eta_all(iM));
    save(fullfile(save_path,save_file),...
         'W','H','objective','iter_times','params');
end


end


