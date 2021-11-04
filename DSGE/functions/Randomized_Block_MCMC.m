function [sampler ln_posterior c acceptance_rate] = Randomized_Block_MCMC(hyperparams,param_name,prior_information,theta_initial,Sigma,Y)

%{
The function draws a sampler from a given posterior
with randomized block MCMC method.
 - inputs
N: number of draws
N_block: number of randomized block
c: scale parameter for Sigma
cb: scale parameter for block adjustment
theta_initial: initial point to start MCMC
Sigma: proposal density covariance
Y: observation

%}

% declare some variables
N = hyperparams.N;
N_block = hyperparams.N_block;
c0 = hyperparams.c0;
target_rate = hyperparams.AccptRate;
constant_c = hyperparams.constant_c;
ME_scale = hyperparams.ME_scale;

% determine burn-in size: at least 1e3, at most 1e4
burnin = max(1e4, min(.20*N, 1e4));

% measure the parameter dimension
n = length(theta_initial);

% estimation basket
sampler = theta_initial;
ln_posterior = zeros(N+burnin+1,1);
c = zeros(N+burnin+1,1);
acceptance_rate = zeros(N+burnin,1);
c(1) = c0;
accepted = 0;

ln_pr0 = log_prior(prior_information,theta_initial);
ln_lkhd0 = log_likelihood_kalman(theta_initial,Y,ME_scale);
ln_pst0 = ln_pr0 + ln_lkhd0;
ln_posterior(1) = ln_pst0;

% clock in
tic;

% MCMC iteration
for i = 1:burnin+N
    
    % generate random block
    [block1 block2] = generate_random_block(n,N_block);
    
    % draw a proposal
    if length(block1) == n
        theta1 = mvnrnd(sampler(:,end),c(i)^2*Sigma)';
    else
        theta1 = block_Markov_chain_proposal(sampler(:,end),Sigma,c(i),block1);
    end
    
    % compute the log posterior of proposal
    ln_pr1 = log_prior(prior_information,theta1);
    ln_lkhd1 = log_likelihood_kalman(theta1,Y,ME_scale);
    ln_pst1 = ln_pr1 + ln_lkhd1;
    
    % acceptance probablity
    acceptance_probability = exp(ln_pst1 - ln_pst0);
    
    % sampler decision
    if acceptance_probability > rand
        % bestow the survivorship
        survivor = theta1;
        % carry the log posterior of winner
        ln_pst0 = ln_pst1;
        % count up the accepted
        accepted = accepted + 1;
    else
        % keep the draw
        survivor = sampler(:,end);
    end
    
    % blockwise iteration
    for b = 2:N_block
        
        % draw a proposal
        thetab = block_Markov_chain_proposal(survivor,Sigma,c(i),block2(b-1,:));
        
        % compute the log posterior of proposal
        ln_pr1 = log_prior(prior_information,thetab);
        ln_lkhd1 = log_likelihood_kalman(thetab,Y,ME_scale);
        ln_pst1 = ln_pr1 + ln_lkhd1;
        
        % acceptance probablity
        acceptance_probability = exp(ln_pst1 - ln_pst0);
        
        % sampler decision
        if acceptance_probability > rand
            % bestow the survivorship
            survivor = thetab;
            % carry the log posterior of winner
            ln_pst0 = ln_pst1;
            % count up the accepted
            accepted = accepted + 1;
        end
        
    end
    
    % compute acceptance rate
    acceptance_rate(i) = accepted / N_block / i;
    
    % adjust scaling factor
    if constant_c == 1
        target_rate = acceptance_rate(i);
    end
    c(i+1) = c(i)*c_adjust(acceptance_rate(i),target_rate);
    
    % add up the survivor to sampler
    sampler = [sampler survivor];
    
    % record log posterior
    ln_posterior(i+1) = ln_pst0;
    
    % interim report
    if mod(i,500) == 0        
        fprintf('RB-MCMC %i/%i drawed(%.2f%s)\n',i,N+burnin,i/(N+burnin)*100,'%')
        fprintf('Acceptance rate: %.2f%s\n',acceptance_rate(i)*100,'%')
        fprintf('---------------------------------\n')
        fprintf('%-15s%5s%11s\n','param','mean','std')
        fprintf('---------------------------------\n')
        fprintf('%-15s%0.4f%12.4f\n',[string(param_name(:,1)),mean(sampler')',std(sampler')']')
        fprintf('_________________________________\n')        
    end

end

sampler = sampler(:,end-N+1:end);

% clock out
toc;

end