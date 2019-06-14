% 2019 04 26  Check param search results


if ismac && isunix  % if on Mac
    addpath /Users/wu-jung/code_git/ooi_sonar/palm_nmf_matlab/
    data_path = '/Users/wu-jung/code_git/ooi_sonar/sample_data/';
    save_path = '/Users/wu-jung/code_git/ooi_sonar/decomp_results/';
elseif isunix  % if on Linux
    addpath ~/code_git/ooi_sonar/palm_nmf_matlab/
    data_path = '~/code_git/ooi_sonar/sample_data/';
    save_path = '~/code_git/ooi_sonar/decomp_results/';
end


rank_all = 1:30;

for rr = 1:length(rank_all)
    fname = sprintf('%s_r%02d_betaW%2.2f_betaH%2.2f_smoothness%08.2f.mat', ...
        'rank_search_20190426', rank_all(rr), ...
        0.1, 0.1, 10);
    load(fullfile(save_path, fname));
    objective_min_10(rr) = objective(end);
    fname = sprintf('%s_r%02d_betaW%2.2f_betaH%2.2f_smoothness%08.2f.mat', ...
        'rank_search_20190426', rank_all(rr), ...
        0.1, 0.1, 100);
    load(fullfile(save_path, fname));
    objective_min_100(rr) = objective(end);
    fname = sprintf('%s_r%02d_betaW%2.2f_betaH%2.2f_smoothness%08.2f.mat', ...
        'rank_search_20190426', rank_all(rr), ...
        0.1, 0.1, 1000);
    load(fullfile(save_path, fname));
    objective_min_1000(rr) = objective(end);
end