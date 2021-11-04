function lprd = log_prior(X,theta)
% compute the log prior density
%Beta(a,b): a = mu^2*(1-mu)/sigma^2 - mu, b = a/mu -a;
%Gamma(a,b): a = (mu/sigma)^2, b = sigma^2 / mu;
%inverse gamma(a,b): a = 2 + (mu/sigma)^2, b = mu*(1+(mu/sigma)^2)
%uniform(a,b): a = mu - sqrt(3)*sigma, b = mu + sqrt(3)*sigma

% measure the parameter dimension
n = length(X);

% start with zero density
lprd = 0;

% add up prior density in turn
for j = 1:n
    mu = X(j,2);
    sigma = X(j,3);
    if X(j,1) == 1 % Beta
        a = mu^2*(1-mu)/sigma^2 - mu;
        b = a/mu -a;
        lprd = lprd + log(betapdf(theta(j),a,b));
    elseif X(j,1) == 2 % Gamma
        a = (mu/sigma)^2;
        b = sigma^2 / mu;
        lprd = lprd + log(gampdf(theta(j),a,b));
    elseif X(j,1) == 3 % inv Gamma
        a = 2 + (mu/sigma)^2;
        b = mu*(1+(mu/sigma)^2);
        lprd = lprd + log(invgampdf(theta(j),a,b));
    elseif X(j,1) == 4 % Normal
        lprd = lprd + log(normpdf(theta(j),mu,sigma));
    elseif X(j,1) == 5 % Uniform
        a = mu - sqrt(3)*sigma;
        b = mu + sqrt(3)*sigma;
        lprd = lprd + log(1/(b-a)*(theta(j)>=a)*(theta(j)<=b));
    end
end
end