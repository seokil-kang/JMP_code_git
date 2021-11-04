load parameter_name
for j = 1:length(nm)    
    eval([char(nm(j)) ' =theta(j);']);
end

% log utility for consumption
chi_c = 1;

% gamma 100 to 1
gamm = gamm100 / 100;

% exponential growth
egamma = exp(gamm);

% time discount rate
Beta = .99;

% adjusted discount rate
betat = Beta * exp(gamm * (1-chi_c));

% leisure preference parameter
vartheta = 1;

% indicator function for internal habit(external if zero)
iota_c = 1;

% steady state depreciation rate
delta = .025;

% production share of capital
Alpha = 1/3;

% steady state elasticity of substitution of differentiated labor
mu_w = .14;

% steady state elasticity of substitution of intermediate goods
mu_p = mu_w;

% average duration of gov't debt
average_duration = 4*5;

% lump-sum transfer is set identified...
rho_z = .98;
varrho_z = .0241;

% steady state tax rates
taunss = .186;
taukss = .218;
taucss = .023;

% steady state target market value of debt share to gdp
s_b = .8;

% steady state gov't spending share to gdp
s_g = .0992;

% measured gross interest rate
R_bar = (egamma / betat - 1) * 100;

% steady state computation
kappa_c = betat / egamma * eta * iota_c;
kappa_i = (1 + betat) * varsigma * egamma^2;
kappa_w = (1 - zeta_w * betat) * (1 - zeta_w) / (1 + betat)...
    / zeta_w / (1 + (1 + mu_w)/mu_w * chi_n);
kappa_p = (1 - zeta_p * betat) * (1 - zeta_p) / (1 + iota_p * betat) / zeta_p;
rho_L = (1 - 1/average_duration) / betat;
rkss = (egamma / betat - 1 + delta) / (1 - taukss);
wss = (1/(1 + mu_p) * (1-Alpha)^(1-Alpha) * Alpha^Alpha *rkss^(-Alpha))^(1/(1-Alpha));
kn = wss / rkss * Alpha / (1 - Alpha);
yn = rkss * kn + wss;
Omegan = kn^Alpha - yn;
ivn = (egamma - 1 + delta) * kn;
cn = (1 - s_g) * yn - ivn;
nss = (1/vartheta * wss / (1 + mu_w) * (1 - taunss) / (1 + taucss) * (1 - kappa_c) / (1 - eta / egamma)*cn^(-chi_c))^(1/(chi_c + chi_n));
kss = kn * nss;
yss = yn * nss;
Omega = Omegan * nss;
ivss = ivn * nss;
css = cn * nss;
tauss = taucss * css + taunss * wss * nss + taukss * rkss * kss;
gss = s_g * yss;
bss = s_b * yss;
zss = (1-1/betat)*bss + tauss - gss;