function ssNMF_runner_normvar_repeat(param_file, row_num, varargin)

% Use parameters defined in a row of a text file to run ss-NMF
%
%   param_file  file storing the input params
%   row_num     row number of param file to be used
%
%   Each row in param_file is a sequence of 6 numbers:
%   rank, betaH, betaW, smoothness, sparsity, max_iter
%
%   varargin{1} is the save_path. If not supplied,
%   save results in the current directory.
%
%   varargin{2} is opt_autostop. If 1-autostop, 
%   if 0-run to max_iter



% Parse input parameters
fileID = fopen(param_file);
P = textscan(fileID,'%d %f %f %f %f %f %d',1,'HeaderLines',row_num);  % read 1 line
fclose(fileID);

r = P{1};
betaH = P{2};
betaW = P{3};
sp = P{4};        % sparsity
sm = P{5};        % smoothness
max_iter = P{6};  % max iteration
rep_num = P{7};      % rep number

rng(rep_num);        % rng seed

if nargin>=3
    save_path = varargin{1};  % path to save results
else
    save_path = pwd;  % if not specified, save results to current path
end

if nargin==4
    opt_autostop = varargin{4};   % is opt_autostop is set
else
    opt_autostop = 0;   % no autostop, run to max_iter
end

rng_start = rng;  % save a copy of the initial rng state


% Set various paths
addpath /home/ch153/wjl/ooi_sonar/palm_nmf_matlab
data_path = '/home/ch153/wjl/ooi_sonar/sample_data';
data_file = '20150817-20151017_MVBS_PCPcleaned.h5';


% If save_path does not exist, create it
if ~exist(save_path, 'dir')
    mkdir(save_path)
end


cd(data_path)  % has to switch folder to ensure reading h5 files correctly


% Get current script name
pname = mfilename('fullpath');
[~,sname,~] = fileparts(pname);  % script name


% Load PCP-cleaned data
L = h5read(data_file,'/L');
LL = L-min(L(:));  % make it non-negative


% Normalize variance across pixels
LL_norm = nan(size(LL));
sigma_all = nan(size(LL,1),1);
for irow = 1:size(LL,1)
    sigma = std(LL(irow,:));
    LL_norm(irow,:) = LL(irow,:)/sigma;
    sigma_all(irow) = sigma;
end
LL = LL_norm;  % use normalized data for decomposition


% Iteration to save varies by order of magnitude
iter_order = 0:4;
iter_save = repmat(1:0.2:9.8,length(iter_order),1);
for iorder = iter_order+1
    iter_save(iorder,:) = iter_save(iorder,:)*10^iter_order(iorder);
end
iter_save = iter_save';
iter_save = unique(int32(iter_save(:)));
iter_save = iter_save(iter_save<=max_iter);
if iter_save(end)~=max_iter
    iter_save(end+1) = max_iter;
end


% Run ss-NMF
params.r = r;
params.betaH = betaH;
params.betaW = betaW;
params.smoothness = sm;
params.sparsity = sp;
params.max_iter = max_iter;
params.opt_autostop = opt_autostop;

fprintf('  rank       = %d\n', r);
fprintf('  betaH      = %.2e\n', betaH);
fprintf('  betaW      = %.2e\n', betaW);
fprintf('  sparsity   = %.2e\n', sp);
fprintf('  smoothness = %.2e\n', sm);
fprintf('  max_iter   = %.2e\n', max_iter);
fprintf('  rep_num    = %d\n', rep_num);
fprintf('  autostop   = %d\n', opt_autostop);

save_file = ...
    sprintf('%s_r%02d_betaW%2.2f_betaH%2.2f_sm%0.2e_sp%0.2e_maxiter%0.2e_rep%03d.mat', ...
            sname, params.r, ...
            params.betaW, params.betaH, ...
            params.smoothness, params.sparsity,max_iter,rep_num);

fprintf('%s\n', save_file);
[W, H, objective, iter_times, iter_save_out, W_init, H_init, W_steps, H_steps] = ...
    palm_nmf_detail(LL, params, iter_save);

save_file = fullfile(save_path, save_file);
m = matfile(save_file,'writable',true);
m.W = W;
m.H = H;
m.W_init = W_init;
m.H_init = H_init;
m.W_steps = W_steps;
m.H_steps = H_steps;
m.objective = objective;
m.iter_save_out = iter_save_out;
m.iter_times = iter_times;
m.iter_save = iter_save;
m.params = params;
m.rng = rng_start;   % starting state of rng
m.sigma_all = sigma_all;

fprintf('%s\n', datetime('now','Format','y-M-d HH:mm:ss'));

