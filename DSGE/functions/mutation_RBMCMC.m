function [particle_mutated ln_posterior ln_likelihood acceptance_rate] = mutation_RBMCMC(hyperparams,prior_information,particle_initial,ln_pst0,ln_lkhd0,mu,Sigma,Y,phi)

%{
The function mutates a particle with randomized block MCMC method.
 - inputs
N_mutate: number of mutation steps
N_block: number of randomized block
c: scale parameter for Sigma
theta_initial: initial point to start MCMC
Sigma: proposal density covariance
Y: observation

%}

% declare some variables
N_mutate = hyperparams.N_mutate;
N_block = hyperparams.N_block;
alpha = hyperparams.alpha;
c = hyperparams.c;
ME_scale = hyperparams.ME_scale;
phi0 = phi(1);
phi1 = phi(2);

% measure the parameter dimension
n = length(particle_initial);

% estimation basket
particle_mutated = particle_initial;
accepted = 0;
%ln_pr0 = log_prior(prior_information,particle_initial);
%ln_lkhd0 = log_likelihood_kalman(particle_initial,Y,ME_scale);
ln_pst0 = ln_pst0 + (phi1 - phi0) * ln_lkhd0;

% MCMC iteration
for i = 1:N_mutate
    
    % generate random block
    [block1 block2] = generate_random_block(n,N_block);
    
    % draw a proposal
    if length(block1) == n
        particle1 = alpha*mvnrnd(particle_mutated,c^2*Sigma)'...
            + (1-alpha)/2 * mvnrnd(particle_mutated,c^2*diag(diag(Sigma)))'...
            + (1-alpha)/2 * mvnrnd(mu,c^2*Sigma)';
    else
        particle1 = alpha*block_Markov_chain_proposal(particle_mutated,Sigma,c,block1)...
            + (1-alpha)/2 * block_Markov_chain_proposal(particle_mutated,diag(diag(Sigma)),c,block1)...
            + (1-alpha)/2 * block_Markov_chain_proposal(mu,Sigma,c,block1);
    end
    
    % compute the log posterior of proposal
    ln_pr1 = log_prior(prior_information,particle1);
    ln_lkhd1 = log_likelihood_kalman(particle1,Y,ME_scale);
    ln_pst1 = ln_pr1 + phi1 * ln_lkhd1;
    
    % acceptance probablity
    acceptance_probability = exp(ln_pst1 - ln_pst0);
    
    % sampler decision
    if acceptance_probability > rand
        % bestow the survivorship
        particle_mutated = particle1;
        % carry the log posterior and likelihood of winner
        ln_pst0 = ln_pst1;
        ln_lkhd0 = ln_lkhd1;
        % count up the accepted
        accepted = accepted + 1;
    end
    
    % blockwise iteration
    for b = 2:N_block
        
        % draw a proposal
        particleb = alpha*block_Markov_chain_proposal(particle_mutated,Sigma,c,block2(b-1,:))...
            +(1-alpha)/2*block_Markov_chain_proposal(particle_mutated,diag(diag(Sigma)),c,block2(b-1,:))...
            +(1-alpha)/2*block_Markov_chain_proposal(mu,Sigma,c,block2(b-1,:));
        
        % compute the log posterior of proposal
        ln_pr1 = log_prior(prior_information,particleb);
        ln_lkhd1 = log_likelihood_kalman(particleb,Y,ME_scale);
        ln_pst1 = ln_pr1 + phi1 * ln_lkhd1;
        
        % acceptance probablity
        acceptance_probability = exp(ln_pst1 - ln_pst0);
        
        % sampler decision
        if acceptance_probability > rand
            % bestow the survivorship
            particle_mutated = particleb;
            % carry the log posterior of winner
            ln_pst0 = ln_pst1;
            ln_lkhd0 = ln_lkhd1;
            % count up the accepted
            accepted = accepted + 1;
        end
        
    end
    
end

% compute the acceptance rate
acceptance_rate = accepted / N_mutate / N_block;

% compute log posterior and likelihood
ln_posterior = ln_pst0;
ln_likelihood = ln_lkhd0;

end