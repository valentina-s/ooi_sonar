% 2019 04 26  Check param search results


if ismac && isunix  % if on Mac
    addpath /Users/wu-jung/code_git/ooi_sonar/
    data_path = '/Users/wu-jung/code_git/ooi_sonar/sample_data/';
    save_path = '/User/wu-jung/code_git/ooi_sonar/decomp_results/';
elseif isunix  % if on Linux
    addpath ~/code_git/ooi_sonar/
    data_path = '~/code_git/ooi_sonar/sample_data/';
    save_path = '~/code_git/ooi_sonar/decomp_results/';
end


rank_all = 1:20;

for rr = 1:length(rank_all)
    fname = sprintf('%s_r%02d_betaW%2.2f_betaH%2.2f_smoothness%6.2f.mat', ...
        'param_search_20190426', rank_all(rr), ...
        0.1, 0.1, 100);
    load(fullfile(save_path, fname));
    objective_min(rr) = objective(end);
end