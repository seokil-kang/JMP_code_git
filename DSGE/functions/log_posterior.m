function ln_pst = log_posterior(theta,Y,prior_info,ME_scale)
% computes log posterior value

% log prior
ln_pr = log_prior(prior_info,theta);
% log likelihood
ln_lkhd = log_likelihood_kalman(theta,Y,ME_scale);

% log posterior
ln_pst = ln_pr + ln_lkhd;

end