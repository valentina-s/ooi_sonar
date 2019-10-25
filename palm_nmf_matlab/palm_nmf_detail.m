% Majority of code lifted from https://github.com/raimon-fa/palm-nmf/blob/master/palm_nmf.m
% Modified so that can save detailed output for individual steps (June 2019)
% Added the 'random_seed' flag since it was not actually implemented before even though option exists (Oct-25 2019)

function [W, H, objective, iter_times, W_init, H_init, W_steps, ...
	   H_steps, rng_save] = palm_nmf_detail(V, params, varargin)
% Algorithm for NMF with eucidian norm as objective function and 
%L1 constraint on W for sparse paterns and Tikhonov regularization 
%for smooth activation coefficients.
%
% [W, H, objective] = sparse_nmf(v, params)
% 
% Objective function: 
% || W * H - V ||² + lambda * || W ||_1 + eta || H T ||² + betaW ||W||²
%                                                        + betaH ||H||²
%
%
% Inputs:
% V:  matrix to be factorized
% params: optional parameters
%     sparsity: weight for the L1 sparsity penalty (default: 0)
%
%     smoothness: weight for the smoothness constraint.
%
%     max_iter: maximum number of iterations (default: 100)
%
%     conv_eps: threshold for early stopping (default: 0, 
%                                             i.e., no early stopping)
%     random_seed: set the random seed, if not set, use rng('default')
%                  values should be allowed inputs of rng
%                  for cluster jobs, set random_seed='shuffle' for true randomization
%                   
%     init_W:   initial setting for W (default: random; 
%                                      either init_w or r have to be set)
%     r: (K)       # basis functions (default: based on init_w's size;
%                                  either init_w or r have to be set)
%     init_H:   initial setting for H (default: random)
%     init_W:   initial setting for W (default: random)
%
%     gamma1:   constant > 1 for the gradient descend step of W.
%
%     gamma2:   constant > 1 for the gradient descend step of W.
%
%     betaH:   constant. L-2 constraint for H.
%     betaW:   constant. L-2 constraint for W.

% Outputs:
% W: matrix of basis functions
% H: matrix of activations
% objective: objective function values throughout the iterations
% iter_times: time passed until iteration ith
% W_init: initial W
% H_init: initial H
% W_steps: W in all iterations
% H_steps: H in all iterations
% rng_save: random number generator state


% Check input for iterations to save
if nargin==3  % if iteration steps are specified
    iter_save = varargin{1};
    fprintf('Iterations to save are: %s\n',...
            strjoin(cellstr(num2str(iter_save(:))),','));
else
    iter_save = [];  % save all steps
    disp('Save all iterations')
end

% Get other parameters
m = size(V, 1);
n = size(V, 2);

if ~exist('params', 'var')
    params = struct;
end
if ~isfield(params, 'max_iter')
    params.max_iter = 100;
end
if ~isfield(params, 'sparsity')
    params.sparsity = 0;
end
if ~isfield(params, 'smoothness')
    params.smoothness = 0;
end
if ~isfield(params, 'conv_eps')
    params.conv_eps = 0;
end
if ~isfield(params, 'betaH')
    params.betaW= 0.1;
end
if ~isfield(params, 'betaH')
    params.betaH= 0.1;
end
if ~isfield(params, 'gamma1')
    gamma1 = 1.001;
else
    gamma1 = params.gamma1;
end
if ~isfield(params, 'gamma2')
    gamma2 = 1.001;
else
    gamma2 = params.gamma2;
end

% Deal with random number generation
if ~isfield(params, 'random_seed')
    params.random_seed = 'default';
end
rng(params.random_seed);
rng_save = rng;

% Initialize W and H
if ~isfield(params, 'init_W')
    if ~isfield(params, 'r')
        error('Number of components or initialization must be given')
    end
    r = params.r;
    W = rand(m, r);
else
    r = size(params.init_W, 2);
    W = params.init_W;
end
if ~isfield(params, 'init_H')
    H = rand(r, n);
else
    H = params.init_H;
end

%%% PALM NMF %%%
lambda = params.sparsity;
eta = params.smoothness;
betaH = params.betaH;
betaW = params.betaW;
objective = zeros(params.max_iter,1);
iter_times = zeros(params.max_iter,1);

%%% Initialize outputs for saving iterations
H_init = H;  % initial H
W_init = W;  % initial W
if ~isempty(iter_save)  % only save specified iterations
    H_steps = zeros(r, n, length(iter_save));  % H size: [m x r]
    W_steps = zeros(m, r, length(iter_save));  % W size: [r x n]
else
    H_steps = zeros(r, n, params.max_iter);  % H size: [m x r]
    W_steps = zeros(m, r, params.max_iter);  % W size: [r x n]
end
save_iter_count = 1;  % counter for saved iterations

tic
if lambda == 0 && eta == 0
    %%% Original NMF %%%
    %In the case without constraints it can be shown that 
    % the gammas can be divided by 2 (Bolte 2014)
    gamma1 = gamma1 / 2;
    gamma2 = gamma2 / 2;
    for it = 1:params.max_iter
        % 1. W updates %%%%%%%%%%%%%%%%%%%%%%%%%%%
        c = gamma1 * 2 * norm(H*H','fro');    
        % gradient descend
        z1 = W - (1 / c) * 2 * ((W * H - V) * H');
        % proximity operator
        z1 = max(z1, 0);
        W = z1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2. H update %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        d = gamma2 * 2 * (norm(W*W','fro'));
        % gradient descend
        z2 = H - (1 / d) * 2 * (W' * (W * H - V));   
        % proximity operator 
        z2 = max(z2,0);
        H = z2;    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute the objective function %%%%%%%%%%
        objective(it) = norm(W * H - V,'fro')^2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        iter_times(it) = toc;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Save H and W step
        if ~isempty(iter_save)
            if ismember(it, iter_save)
                H_steps(:,:,save_iter_count) = H;
                W_steps(:,:,save_iter_count) = W;
                save_iter_count = save_iter_count+1;
            end
        else
            H_steps(:,:,it) = H;
            W_steps(:,:,it) = W;
        end
    end
elseif lambda > 0 && eta == 0
    %%% SPARSE NMF %%%
    for it = 1:params.max_iter
        % 1. W updates %%%%%%%%%%%%%%%%%%%%%%%%%%%
        c = gamma1 * 2 * norm(H*H','fro');    
        % gradient descend
        z1 = W - (1 / c) * 2 * ((W * H - V) * H');
        % proximity operator
        z1 = max(z1 - 2 * lambda / c, 0);
        W = z1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2. H update %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        d = gamma2 * 2 * (norm(W*W','fro') + betaH);
        % gradient descend
        z2 = H - (1 / d) * 2 * (W' * (W * H - V) + betaH * H);   
        % proximity operator 
        z2(z2 < 0) = 0;
        H = z2;    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute the objective function %%%%%%%%%%
        objective(it) = norm(W * H - V,'fro')^2 + ...
            lambda * sum(sum(abs(W))) + betaH * norm(H,'fro');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        iter_times(it) = toc;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Save H and W step
        if ~isempty(iter_save)
            if ismember(it, iter_save)
                H_steps(:,:,save_iter_count) = H;
                W_steps(:,:,save_iter_count) = W;
                save_iter_count = save_iter_count+1;
            end
        else
            H_steps(:,:,it) = H;
            W_steps(:,:,it) = W;
        end
    end
elseif lambda == 0 && eta > 0
    %%% SMOOTH NMF %%%
    %Tikhonov regularization matrix
    T = eye(n) - diag(ones(n-1,1),-1);
    T = T(:,1:end-1);
    TTp = T*T';
    TTp_norm = norm(TTp,'fro');
    for it = 1:params.max_iter
        % 1. W updates %%%%%%%%%%%%%%%%%%%%%%%%%%%
        c = gamma1 * 2 * (norm(H*H','fro') + betaW);    
        % gradient descend
        z1 = W - (1 / c) * 2 * ((W * H - V) * H' + betaW * W);
        % proximity operator 
        z1 = max(z1, 0);
        W = z1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2. H update %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        d = gamma2 * 2 * (norm(W*W','fro') + eta * TTp_norm);
        % gradient descend
        z2 = H - (1 / d) * 2 * (W' * (W * H - V) + eta * (H * TTp));   
        % proximity operator 
        z2 = max(z2,0);
        H = z2;    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute the objective function %%%%%%%%%%
        objective(it) = norm(W * H - V,'fro')^2 + ...
            eta * norm(H * T,'fro')^2 + betaW * norm(W,'fro'); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        iter_times(it) = toc;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Save H and W step
        if ~isempty(iter_save)
            if ismember(it, iter_save)
                H_steps(:,:,save_iter_count) = H;
                W_steps(:,:,save_iter_count) = W;
                save_iter_count = save_iter_count+1;
            end
        else
            H_steps(:,:,it) = H;
            W_steps(:,:,it) = W;
        end
    end  
elseif lambda > 0 && eta > 0
    %%% SMOOTH and SPARSE NMF %%%
    %Tikhonov regularization matrix
    T = eye(n) - diag(ones(n-1,1),-1);
    T = T(:,1:end-1);
    TTp = T*T';
    TTp_norm = norm(TTp,'fro');
    for it = 1:params.max_iter
        % 1. W updates %%%%%%%%%%%%%%%%%%%%%%%%%%%
        c = gamma1 * 2 * (norm(H*H','fro') + betaW);    
        % gradient descend
        z1 = W - (1 / c) * 2 * ((W * H - V) * H' + betaW * W);
        % proximity operator 
        z1 = max(z1 - 2 * lambda / c, 0);
        W = z1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2. H update %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        d = gamma2 * 2 * (norm(W*W','fro') + eta * TTp_norm + betaH);
        % gradient descend
        z2 = H - (1 / d) * 2 * (W' * (W * H - V) + eta * (H * TTp) +...
                                betaH * H);   
        % proximity operator 
        z2 = max(z2,0);
        H = z2;    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute the objective function %%%%%%%%%%
        objective(it) = norm(W * H - V,'fro')^2 + ...
            eta * norm(H * T,'fro')^2 + lambda * sum(sum(abs(W))) + ...
            betaW * norm(W,'fro') + betaH * norm(H,'fro'); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        iter_times(it) = toc;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Save H and W step
        if ~isempty(iter_save)
            if ismember(it, iter_save)
                H_steps(:,:,save_iter_count) = H;
                W_steps(:,:,save_iter_count) = W;
                save_iter_count = save_iter_count+1;
            end
        else
            H_steps(:,:,it) = H;
            W_steps(:,:,it) = W;
        end
    end
else
    error('Give positive values to the parameters')
end
extime = toc;
objective = objective(1:it);
disp(['Number of iterations = ' num2str(it)])
disp(['Execution time = ' num2str(extime) ' sec' ])
end
