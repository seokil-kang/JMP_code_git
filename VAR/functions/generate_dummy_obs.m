function [Yd Xd] = generate_dummy_obs(lambda,Y_full,X,S_prior,A_prior,c_prior)

% measure the data and model dimension
[T K] = size(X);
[T_full, M] = size(Y_full);
p = T_full - T;

% construct the dependent variable
Y = Y_full(p+1:end,:);

% some initial measurements
y0bar = mean(Y_full(1:p,:))';

% 1st block: coefficients
Yd = [lambda(1) * A_prior * S_prior; zeros(M*(p-1),M)];
Xd = [zeros(M*p,1) kron(lambda(1)*diag(1:p)^lambda(2),S_prior)];

% 2nd block: constant terms
Yd = [Yd; lambda(3)*c_prior'];
Xd = [Xd; lambda(3) zeros(1,M*p)];

% 3rd block: co-persistence(initial value conditionality)
Yd = [Yd; lambda(4)*y0bar'];
Xd = [Xd; lambda(4) kron(ones(1,p),lambda(4)*y0bar')];

% 4th block: sum-of-coefficients(no cointegration)
Yd = [Yd; lambda(5)*diag(y0bar)];
Xd = [Xd; zeros(M,1) kron(ones(1,p),lambda(5)*diag(y0bar))];
    
% 5th block: covariance mean for diagonal terms
Yd = [Yd; kron(ones(lambda(6),1),S_prior)];
Xd = [Xd; zeros(M*lambda(6),M*p+1)];

end