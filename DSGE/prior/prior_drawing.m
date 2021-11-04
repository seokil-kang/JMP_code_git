function theta = prior_drawing(X,N)
% random draw parameters based the prior_info X with N times
%Beta(a,b): a = mu^2*(1-mu)/sigma^2 - mu, b = a/mu -a;
%Gamma(a,b): a = (mu/sigma)^2, b = sigma^2 / mu;
%inverse gamma(a,b): a = 2 + (mu/sigma)^2, b = mu*(1+(mu/sigma)^2)
%uniform(a,b): a = mu - sqrt(3)*sigma, b = mu + sqrt(3)*sigma

% measure the parameter dimension
[n,~] = size(X);


% draw parameters in turn
theta = [];
for j = 1:n
    mu = X(j,2);
    sigma = X(j,3);
    if X(j,1) == 1 % Beta
        a = mu^2*(1-mu)/sigma^2 - mu;
        b = a/mu -a;
        theta = [theta betarnd(a,b,[N 1])];
    elseif X(j,1) == 2 % Gamma
        a = (mu/sigma)^2;
        b = sigma^2 / mu;
        theta = [theta gamrnd(a,b,[N 1])];
    elseif X(j,1) == 3 % inv Gamma
        a = 2 + (mu/sigma)^2;
        b = mu*(1+(mu/sigma)^2);
        theta = [theta invgamrnd(a,b,[N 1])];
    elseif X(j,1) == 4 % Normal
        theta = [theta normrnd(mu,sigma,[N 1])];
    elseif X(j,1) == 5 % Uniform
        a = mu - sqrt(3)*sigma;
        b = mu + sqrt(3)*sigma;
        theta = [theta rand(N,1)*(b-a) + a];
    end
end
theta = theta';
end