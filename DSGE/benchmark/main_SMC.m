%Sequential Monte Carlo
%Seokil Kang(sk86@iu.edu), Nov 2020

%--------------------
% pre-estimation part
%--------------------

% house keeping
clear; clc; close all;
delete *.asv

%-------------------------------------------------------
% estimation parameter selection
N = 1e4; % # of particles
N_stage = 500; % # of SMC stages
N_mutate = 1; % mutation steps
N_block = 6; % # of block for mutation
alpha = .9; % mixture ratio for proposal density
c0 = 1; % initial scaling factor
target_AccptRate = .25; % target acceptance rate for mutation
ME_scale = 1; % measurement error ratio
lambda = 2.1; % tempering schedule parameter
N_threshold = N/2; % resampling threshold
strata = 1; % resampling method: 1(strata) 0(systematic)
%-------------------------------------------------------

% keep the default path
path_default = path;

% add function path
addpath('../../data','../functions','../gensys','model','../prior','../sims_minwel')

% call prior parameter
[prior_info param_name] = prior_info_setup;

% measure the parameter dimension
n = length(prior_info);

% load observation
% data_index: vector of integers with the following order
% 1 2 3 4 5 6 7 8 9 10 11 12 13 14
% y c i w n b s g z p  R  Rb x s(accounting)
[y dateq data_index diff_index] = data_importing([1 3 4 5 6 7 10 11 12],'1960-1','2006-4', 1);

% estimation baskets
theta = zeros(n,N,N_stage+1); % param particle with stage 0
theta_hat = zeros(n,N,N_stage); % selected param particle
mutated = zeros(n,N); % mutation basket for parallel for loop
W = zeros(N,N_stage+1); % particle weights with stage 0
w_tilde = zeros(N,N_stage); % incremental weights
ln_likelihood = zeros(N,1); % likelihood basket
ln_posterior = zeros(N,N_stage+1); % log posterior for each particle
rho = zeros(N_stage,1); % resampling indicator
ESS = zeros(N_stage,1); % effective sample size
Acceptance_rate = zeros(N_stage,1); % acceptance rate
c = zeros(N_stage+1,1); % scaling factor
c(1) = c0;

% construct the hyper parameter structure as a single input for mutation
hyperparams_mutation = struct('N_mutate',N_mutate,'N_block',N_block,'alpha',alpha,'c',c(1),'ME_scale',ME_scale);

% tempering schedule
phi = linspace(0,1,N_stage+1)'.^lambda;

% SMC algorithm
fprintf('Sequential Monte Carlo begins...\n')

%% ----------------
% 1. Initialization
%------------------

% clock in
tic;

% turn on the parallel computing mode
myCluster = parcluster('local')
parpool(myCluster.NumWorkers)

fprintf('_________________________________\n')
fprintf('SMC prior stage of %i particles\n',N)

% Draw the initial particles from the prior
theta(:,:,1) = prior_drawing(prior_info,N);

% compute the likelihood for the Correction step
parfor j = 1:N
    ln_likelihood(j) = log_likelihood_kalman(theta(:,j,1),y,ME_scale);
    ln_posterior(j,1) = log_prior(prior_info,theta(:,j,1)) + ln_likelihood(j);
end

% set the uniform initial weights
W(:,1) = 1;

% compute 1st and 2nd moment of prior
mnt_1st = theta(:,:,1) * W(:,1)/N;
mnt_2nd = (theta(:,:,1) - mnt_1st).^2 * W(:,1)/N;
fprintf('---------------------------------\n')
fprintf('%-15s%5s%11s\n','param','mean','std')
fprintf('---------------------------------\n')
fprintf('%-15s%0.4f%12.4f\n',[string(param_name(:,1)),mnt_1st,sqrt(mnt_2nd)]')
fprintf('_________________________________\n')
%%
%-------------
% 2. Recursion
%-------------

for stage = 1:N_stage
    %---------------
    % (a) Correction
    %---------------
    
    % reweight the particles from the previous stage
    % by defining log of the incremental weights
    ln_w_tilde = (phi(stage+1) - phi(stage)) * ln_likelihood;
    
    % transform the log into the original weight
    % scale issue: log value is too small so exp(ln_w) = 0 usually...
    w_tilde(:,stage) = exp(ln_w_tilde);
    
    % normalize weights
    W_tilde = w_tilde(:,stage) .* W(:,stage) / mean(w_tilde(:,stage) .* W(:,stage));
    
    %--------------
    % (b) Selection
    %--------------
    
    % compute effective sample size
    ESS(stage) = N / (1/N * W_tilde' * W_tilde);
    
    % compute resampling indicator
    rho(stage) = (ESS(stage) < N_threshold);
    
    % crossroads for resampling
    if rho(stage) % resample the particles       
        
        if strata % resampling via stratified sampling method
            xi = rand(N,1);
        else % resampling via systematic method
            xi = rand;
        end
        % create sampling stratum
        u = linspace(0,N-1,N)' + xi;
        
        % local bootstrap resampling
        for i = 1:N
            A(i) = max(1,find(cumsum(W_tilde)>=u(i),1));
        end
        
        % resampled particle
        theta_hat(:,:,stage) = theta(:,A,stage);
        
        % normalize the weight uniformly
        W(:,stage+1) = 1;
                
    else % do not resample the particles        
        theta_hat(:,:,stage) = theta(:,:,stage);
        W(:,stage+1) = W_tilde;
    end
    
    % compute the first moment of particle
    mnt_1st = theta_hat(:,:,stage) * W(:,stage+1) / N;
    
    % compute the second moment of particle
    % two stage because matlab easily breaks symmetry of matrix
    mntdev = (theta_hat(:,:,stage)-mnt_1st) * diag(sqrt(W(:,stage+1)));
    mnt_2nd = mntdev * mntdev'/N;
    
    %-------------
    % (c) Mutation
    %-------------
    
    % run the randomized block MCMC algorithm
    parfor i = 1:N
        [mutatedi ln_psti ln_likehd acpt_rate] = mutation_RBMCMC(...
            hyperparams_mutation,prior_info,...
            theta_hat(:,i,stage),ln_posterior(i,stage),ln_likelihood(i),...
            mnt_1st,mnt_2nd,y,phi(stage:stage+1));
        mutated(:,i) = mutatedi;
        ln_pst(i) = ln_psti;
        ln_likelihood(i) = ln_likehd;
        accepted_stage(i) = acpt_rate;
    end
    % update the particle
    theta(:,:,stage+1) = mutated;
    ln_posterior(:,stage+1) = ln_pst;
    
    % compute the average acceptance rate
    Acceptance_rate(stage) = mean(accepted_stage);
    
    % update the scaling factor
    c(stage+1) = c(stage) * c_adjust(Acceptance_rate(stage),target_AccptRate);
    hyperparams_mutation.c = c(stage+1);
    
    % compute the first moment of particle(after mutation)
    mnt_1st = theta_hat(:,:,stage) * W(:,stage+1) / N;
    
    % compute the second moment of particle(after mutation)
    % two stage because matlab easily breaks symmetry of matrix
    mntdev = (theta_hat(:,:,stage)-mnt_1st) * diag(sqrt(W(:,stage+1)));
    mnt_2nd = mntdev * mntdev'/N;
    
    % stage report
    fprintf('SMC stage %i/%i of %i particles \n',stage,N_stage,N)
    if rho(stage) == 1
        fprintf('particles are resampled!\n')        
    else
        fprintf('particles are not resampled\n')
    end
    fprintf('%i total resampling, ESS = %.2f\n',sum(rho),ESS(stage))
    fprintf('Average mutation = %.2f%s\n',Acceptance_rate(stage)*100,'%')
    fprintf('updated scaling factor = %.2f\n',c(stage+1))
    fprintf('---------------------------------\n')    
    fprintf('%-15s%5s%11s\n','param','mean','std')
    fprintf('---------------------------------\n')
    fprintf('%-15s%0.4f%12.4f\n',[string(param_name(:,1)),mnt_1st,sqrt(diag(mnt_2nd))]')    
    fprintf('_________________________________\n')
    
end

fprintf('SMC completed.\n')

% clock out
toc;

% retrieve path
path(path_default);

% save the dataset
y_data = y;

% save results
save('result/SMC_posterior.mat','theta','W','w_tilde','c','ESS','ln_posterior','Acceptance_rate','y_data','dateq','-v7.3')