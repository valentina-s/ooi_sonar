% 2019 04 26  Check param search results


if ismac && isunix  % if on Mac
    addpath /Users/wu-jung/code_git/ooi_sonar/palm_nmf_matlab/
    data_path = '/Users/wu-jung/code_git/ooi_sonar/sample_data/';
    save_path = '/Users/wu-jung/code_git/ooi_sonar/decomp_results_smoothness/';
elseif isunix  % if on Linux
    addpath ~/code_git/ooi_sonar/palm_nmf_matlab/
    data_path = '~/code_git/ooi_sonar/sample_data/';
    save_path = '~/code_git/ooi_sonar/decomp_results_smoothness/';
end


sm_all = [0,1,2,5,10,20,50,100,200,500,1000,2000,5000,10000];

objective_min = zeros(length(sm_all),1);
for rr = 1:length(sm_all)
    fname = sprintf('%s_r%02d_betaW%2.2f_betaH%2.2f_smoothness%09.2f.mat', ...
        'smoothness_search_20190426', 3, 0.1, 0.1, sm_all(rr));
    load(fullfile(save_path, fname));
    objective_min(rr) = objective(end);
end

figure
semilogx(sm_all, objective_min,'o-','linewidth',2)
set(gca,'fontsize',14)
xlabel('Smoothness');
ylabel('Objective');
title('Influence of smoothness constraints','fontsize',16);
