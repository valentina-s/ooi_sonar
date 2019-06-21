% 2019 06 21  Subsample previously saved iterations to save space

addpath ~/code_git/ooi_sonar/palm_nmf_matlab
%save_path = ['/media/wu-jung/internal_2tb/nmf_results/' ...
%             'decomp_lambda_rank_sweep_revisit_1e4'];
%save_path = ['/media/wu-jung/internal_2tb/nmf_results/' ...
%             'decomp_smoothness_revisit_2e4'];
save_path = ['/media/wu-jung/internal_2tb/nmf_results/' ...
             'decomp_beta_revisit_1e4'];

% ss-NMF params
rank = 3;
sm = [1e2,1e4,1e6,1e8];
beta = 0.0;
sp = [2,5,10,20];
%fname_pre = 'lambda_20190617_rank';
%fname_pre = 'smoothness_search_20190611';
fname_pre = 'beta_20190616';

% Iteration to save varies by order of magnitude
iter_order = 0:4;
iter_save = repmat(1:0.2:9,length(iter_order),1);
for iorder = iter_order+1;
    iter_save(iorder,:) = iter_save(iorder,:)*10^iter_order(iorder);
end
iter_save = iter_save';
iter_save = unique(int32(iter_save(:)));

% Loop through all files
for iloop = 1:length(sm)

    % Filename stuff
    %fname_ori = sprintf('%s_r%02d_betaW%2.2f_betaH%2.2f_smoothness%0.2e_sparsity%0.2e',...
    %                        fname_pre, rank, ...
    %                        beta, beta, ...
    %                        sm, sp(iloop));
    fname_pre_long = sprintf('%s_r%02d_betaW%2.2f_betaH%2.2f_smoothness%09.2f',...
                        fname_pre, rank, ...
                        beta, beta, ...
                        sm(iloop));

    fname_ori = [fname_pre_long,'.mat'];
    fname_old = [fname_pre_long,'_old.mat'];
    
    fprintf('subsampling: %s\n',fname_ori);

    % Changed old filename and load data
    movefile(fullfile(save_path,fname_ori), ...
             fullfile(save_path,fname_old));
    m_old = load(fullfile(save_path, fname_old));
    
    % Adjust subsampling index
    max_iter = length(m_old.objective);
    iter_save = iter_save(iter_save<=max_iter);
    if iter_save(end)~=max_iter
        iter_save(end+1) = max_iter;
    end

    % Subsample and save to original filename
    m = matfile(fullfile(save_path,fname_ori),'Writable',true);
    m.W = m_old.W;
    m.H = m_old.H;
    m.W_init = m_old.W_init;
    m.H_init = m_old.H_init;
    m.W_steps = m_old.W_steps(:,:,iter_save);
    m.H_steps = m_old.H_steps(:,:,iter_save);
    m.objective = m_old.objective;
    m.iter_times = m_old.iter_times;
    m.iter_save = iter_save;  % this is new from subsampling
    m.params = m_old.params;
end