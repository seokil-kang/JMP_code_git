function [D H W V] = measurement_equation(theta)

% construct measurement equation of the state space representation
% input: theta(model parameter)
% output: the function returns the linear system of measurement equation
% Y(t) = D + H * X(t) + W * u(t), u ~ N(0,V)
% W = I by default
% V must be computed outside of the function

model_setup
variable_listing

for measurement_equation_list = 1
    obs_list = {
        'GDP' % Gross Domestic Product
        'FPI' % Fixed private investment
        'CPHNBS' % Compensation per hour of nonfarm business sector
        'NNBS' % Average weekly hours of nonfarm business sector times employment level
        'MVDP' % Market value of debt held by private from Hall et al.(2018)
        'PSGBC' % Primary surplus from government budget constraint
        'GDPDEFLATOR' % GDP implicit price deflator
        'EFFR' % Effective federal funds rate
        'NHPR' % Nominal holding period return of government debt portfolio from Hall et al.(2018)
        };
    n_obs = length(obs_list);
    for obs_index = 1:n_obs
        eval([char(obs_list(obs_index)) ' = obs_index;']);
    end
end

% measurement equation matrix baskets
D = zeros(n_obs,1);
H = zeros(n_obs,n_var);
W = zeros(n_obs,1);
W(PSGBC) = 1;
V = sigma_ME_s;

% output
D(GDP) = gamm100;
H(GDP,y) = 1e2;
H(GDP,Ly) = -1e2;
H(GDP,u_gamma) = 1e2;

% investment
D(FPI) = gamm100;
H(FPI,ivst) = 1e2;
H(FPI,Livst) = -1e2;
H(FPI,u_gamma) = 1e2;

% wage
D(CPHNBS) = gamm100;
H(CPHNBS,w) = 1e2;
H(CPHNBS,Lw) = -1e2;
H(CPHNBS,u_gamma) = 1e2;

% labor
D(NNBS) = n_bar;
H(NNBS,n) = 1e2;

% debt
D(MVDP) = gamm100;
H(MVDP,b) = 1e2;
H(MVDP,Lb) = -1e2;
H(MVDP,u_gamma) = 1e2;

% primary surplus
D(PSGBC) = s_bar - gamm100 + R_bar;
H(PSGBC,s) = 1e2;
H(PSGBC,b) = 1e2*(betat-1);

% inflation
D(GDPDEFLATOR) = pi_bar;
H(GDPDEFLATOR,infl) = 1e2;

% interest rate
D(EFFR) = pi_bar + R_bar;
H(EFFR,R) = 1e2;

% bond return
D(NHPR) = pi_bar + R_bar;
H(NHPR,RL) = 1e2;
end