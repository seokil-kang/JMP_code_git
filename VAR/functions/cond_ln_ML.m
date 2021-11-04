function lnML = cond_ln_ML(lambda,lambda6,Y_full,X,S_prior,A_prior,c_prior)

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
Yd = [Yd; kron(ones(lambda6,1),S_prior)];
Xd = [Xd; zeros(M*lambda6,M*p+1)];

% compute marginal likelihood
Y_star = [Yd; Y];
X_star = [Xd; X];
%A_star = inv(X_star' * X_star) * X_star' * Y_star;
A_star = (Y_star' / X_star')';
S_star = (Y_star - X_star * A_star)' * (Y_star - X_star * A_star);
[T_star,~] = size(Y_star);

Ad = (Yd' / Xd')';
Sd = (Yd - Xd * Ad)' * (Yd - Xd * Ad);
[Td,~] = size(Yd);

%det(X_star'*X_star)^(-M/2) * det(S_star)^(-(T_star-K)/2)*prod(gamma((T_star-K+1-(1:M))/2))
numer = -M/2 * log(det(X_star'*X_star)) -(T_star-K)/2 * log(det(S_star)) + sum(gammaln((T_star-K+1-(1:M))/2));
denom = -M/2 * log(det(Xd'*Xd)) -(Td-K)/2 * log(det(Sd)) + sum(gammaln((Td-K+1-(1:M))/2));
lnML = -M*T/2*log(2*pi) + numer - denom;

end