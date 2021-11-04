function theta_pre_condition = preliminary_condition(theta)
%indicator function for model specific preliminary (boundary) condition
%parameter, in order to save time for estimation. 
%e.g. if a MCMC proposal gives impossible parameter, 
%by concluding its likelihood to zero directly allows one 
%to skip the solution&likeilhood evaluation step

load prior_spec

% parameters whose prior has a positive support
theta_positive = theta(find(prior_spec <= 3));
positive_bound = prod(theta_positive>0);

% parameters whose prior follows Beta bounded below one
theta_beta = theta(find(prior_spec == 1));
beta_bound = prod(theta_beta<1);

model_setup

% fixed cost for production should not be negative
ss_check = prod(Omega >= 0);

% final decision rule
theta_pre_condition = positive_bound * beta_bound * ss_check;
end