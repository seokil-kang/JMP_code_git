function [G0 G1 C Psi Pi S] = equilibrium_condition(theta)
%represent the model equilibrium condition system
%input: theta(model parameter)
%output: system of model fits to gensys

model_setup
variable_listing
equilibrium_conditions_listing

% construct system matrix baskets
% G0*y(t+1) =  G1*y(t) + Psi*shock(t+1) + Pi*forecast_error(t+1) + C
G0 = zeros(n_eqbm_cond,n_var);
G1 = zeros(n_eqbm_cond,n_var);
C = zeros(n_eqbm_cond,1);
Psi = zeros(n_eqbm_cond,n_shock);
Pi = zeros(n_eqbm_cond,n_forecast_error);

% model equilibrium conditions

%-------------------------------------------
% marginal utility of consumption proxy(MUc)
%-------------------------------------------
G0(MUc,Uc) = 1/chi_c;
G0(MUc,u_U) = -1/chi_c;
G0(MUc,c) = egamma / (egamma - eta);
G0(MUc,u_gamma) = eta / (egamma - eta);
G1(MUc,c) = eta / (egamma - eta);

%----------------------------------------
% household FOC w.r.t consumption(HHFOCc)
%----------------------------------------
G0(HHFOCc,Uc) = 1 / (1 - kappa_c);
G0(HHFOCc,xu_gamma) = kappa_c;
G0(HHFOCc,xUc) = -kappa_c;
G0(HHFOCc,lambda) = -1;

%---------------------------------------
% household FOC w.r.t investment(HHFOCi)
%---------------------------------------
G0(HHFOCi,ivst) = 1;
G0(HHFOCi,u_gamma) = 1 / (1 + betat);
G0(HHFOCi,q) = -1 / kappa_i;
G0(HHFOCi,xiv) = -betat / (1 + betat);
G0(HHFOCi,xu_gamma) = -betat / (1 + betat);
G0(HHFOCi,u_i) = -1;
G1(HHFOCi,ivst) = 1 / (1 + betat);

%------------------------------------------------
% household FOC w.r.t capital utilization(HHFOCu)
%------------------------------------------------
G0(HHFOCu,rk) = 1;
G0(HHFOCu,v) = -nu / (1 - nu);

%------------------------------------
% household FOC w.r.t capital(HHFOCk)
%------------------------------------
G0(HHFOCk,q) = -1;
G0(HHFOCk,xlambda) = 1;
G0(HHFOCk,lambda) = -1;
G0(HHFOCk,xu_gamma) = -1;
G0(HHFOCk,xrk) = betat / egamma * (egamma / betat - 1 + delta);
G0(HHFOCk,xq) = betat / egamma * (1 - delta);

%------------------------------------
% intertemporal Euler equation(Euler)
%------------------------------------
G0(Euler,xlambda) = 1;
G0(Euler,xu_gamma) = -1;
G0(Euler,xinfl) = -1;
G0(Euler,lambda) = -1;
G0(Euler,R) = 1;

%--------------------------------------
% longterm govt bond return(BondReturn)
%--------------------------------------
G0(BondReturn,RL) = 1;
G0(BondReturn,Q) = -rho_L * betat / egamma;
G1(BondReturn,Q) = -1;

%----------------------------
% No arbitrage condition(NAC)
%----------------------------
G0(NAC,R) = 1;
G0(NAC,xRL) = -1;

%---------------------------------
% wage markup equation(WageMarkup)
%---------------------------------
G0(WageMarkup,wmk) = -1;
G0(WageMarkup,w) = 1;
G0(WageMarkup,n) = -chi_n;
G0(WageMarkup,u_U) = -1;
G0(WageMarkup,lambda) = 1;

%-------------------------------
% wage Phillips curve(Phillipsw)
%-------------------------------
G1(Phillipsw,w) = 1 / (1 + betat);
G1(Phillipsw,u_gamma) = iota_w / (1 + betat);
G1(Phillipsw,infl) = iota_w / (1 + betat);
G0(Phillipsw,w) = 1;
G0(Phillipsw,xw) = -betat / (1 + betat);
G0(Phillipsw,xinfl) = -betat / (1 + betat);
G0(Phillipsw,xu_gamma) = -betat / (1 + betat);
G0(Phillipsw,infl) = (1 + betat * iota_w) / (1 + betat);
G0(Phillipsw,u_gamma) = (1 + betat * iota_w) / (1 + betat);
G0(Phillipsw,wmk) = kappa_w;
G0(Phillipsw,u_w) = -1;

%--------------------------------
% Law of motion for capital(LOMk)
%--------------------------------
G0(LOMk,ks) = -egamma;
G1(LOMk,ks) = -(1 - delta);
G0(LOMk,u_gamma) = -(1 - delta);
G0(LOMk,ivst) = egamma - 1 + delta;
G0(LOMk,u_i) = (egamma - 1 + delta) * kappa_i;

%------------------------------------
% capital utilization equation(Kutil)
%------------------------------------
G0(Kutil,k) = 1;
G0(Kutil,u_gamma) = 1;
G0(Kutil,v) = -1;
G1(Kutil,ks) = 1;

%-----------------------------------------------
% marginal rate of technology substitution(MRTS)
%-----------------------------------------------
G0(MRTS,w) = 1;
G0(MRTS,n) = 1;
G0(MRTS,rk) = -1;
G0(MRTS,k) = -1;

%-----------------------------------
% price markup equation(PriceMarkup)
%-----------------------------------
G0(PriceMarkup,pmk) = -1;
G0(PriceMarkup,rk) = Alpha;
G0(PriceMarkup,w) = 1 - Alpha;

%--------------------------------
% price Phillips curve(Phillipsp)
%--------------------------------
G0(Phillipsp,infl) = -1;
G0(Phillipsp,xinfl) = betat / (1 + betat * iota_p);
G1(Phillipsp,infl) = -iota_p / (1 + betat * iota_p);
G0(Phillipsp,pmk) = kappa_p;
G0(Phillipsp,u_p) = 1;

%-----------------------------------------------
% aggregate market clearing condition(AggMktClr)
%-----------------------------------------------
G0(AggMktClr,y) = -yss;
G0(AggMktClr,c) = css;
G0(AggMktClr,ivst) = ivss;
G0(AggMktClr,g) = gss;
G0(AggMktClr,v) = (egamma / betat - 1 + delta) * kss;

%------------------------------------------
% aggregate production function(AggProdFcn)
%------------------------------------------
G0(AggProdFcn,y) = -yss / (yss + Omega);
G0(AggProdFcn,k) = Alpha;
G0(AggProdFcn,n) = 1 - Alpha;

%-------------------------
% monetary policy rule(MP)
%-------------------------
G0(MP,R) = -1;
G1(MP,R) = -rho_R;
G0(MP,infl) = (1 - rho_R) * phi_pi;
G0(MP,y) = (1 - rho_R) * phi_x;
G0(MP,u_M) = 1;

%--------------------------
% total tax revenue(TAXREV)
%--------------------------
G0(TAXREV,tau) = -tauss;
G0(TAXREV,c) = taucss * css;
G0(TAXREV,w) = taunss * wss * nss;
G0(TAXREV,n) = taunss * wss * nss;
G0(TAXREV,rk) = taukss * rkss * kss;
G0(TAXREV,k) = taukss * rkss * kss;

%----------------------------------
% primary surplus definition(PSDEF)
%----------------------------------
G0(PSDEF,s) = -bss / betat;
G0(PSDEF,tau) = tauss;
G0(PSDEF,g) = -gss;
G0(PSDEF,z) = -zss;

%----------------------------------
% government budget constraint(GBC)
%----------------------------------
G0(GBC,b) = -betat;
G0(GBC,s) = -1;
G0(GBC,RL) = 1;
G0(GBC,u_gamma) = -1;
G0(GBC,infl) = -1;
G1(GBC,b) = -1;

%--------------------------------------------
% fiscal policy rule: goverment spending(FPg)
%--------------------------------------------
G0(FPg,g) = 1;
G1(FPg,g) = rho_g;
G0(FPg,y) = -(1 - rho_g) * varphi_x;
G1(FPg,b) = -(1 - rho_g) * varphi_b;
G0(FPg,u_g) = -1;

%-------------------------------------------
% fiscal policy rule: lump-sum transfer(FPz)
%-------------------------------------------
G0(FPz,z) = 1;
G1(FPz,z) = rho_z;
G0(FPz,y) = -(1 - rho_z) * psi_x;
G1(FPz,b) = -(1 - rho_z) * psi_b;
G0(FPz,u_z) = -1;

%---------------------------------------
% market value of outstanding debt(MKVB)
%---------------------------------------
G0(MKVB,mktb) = 1;
G0(MKVB,RL) = -1;
G0(MKVB,infl) = 1;
G0(MKVB,u_gamma) = 1;
G1(MKVB,b) = 1;

%------------------------------------------------
% exogenous process for permenant growth(ARgamma)
%------------------------------------------------
G0(ARgamma,u_gamma) = 1;
G1(ARgamma,u_gamma) = varrho_gamma;
Psi(ARgamma,eps_u_gamma) = 1;

%------------------------------------
% exogenous process for int-pref(ARU)
%------------------------------------
G0(ARU,u_U) = 1;
G1(ARU,u_U) = varrho_U;
Psi(ARU,eps_u_U) = 1;

%-------------------------------------------------------------
% exogenous process for marginal efficiency of investment(ARi)
%-------------------------------------------------------------
G0(ARi,u_i) = 1;
G1(ARi,u_i) = varrho_i;
Psi(ARi,eps_u_i) = 1;

%---------------------------------------
% exogenous process for wage markup(ARw)
%---------------------------------------
G0(ARw,u_w) = 1;
G1(ARw,u_w) = varrho_w;
Psi(ARw,eps_u_w) = 1;

%----------------------------------------
% exogenous process for price markup(ARp)
%----------------------------------------
G0(ARp,u_p) = 1;
G1(ARp,u_p) = varrho_p;
Psi(ARp,eps_u_p) = 1;

%-----------------------------------------------
% exogenous process for government spending(ARg)
%-----------------------------------------------
G0(ARg,u_g) = 1;
G1(ARg,u_g) = varrho_g;
Psi(ARg,eps_u_g) = 1;

%-------------------------------------
% exogenous process for labor tax(ARz)
%-------------------------------------
G0(ARz,u_z) = 1;
G1(ARz,u_z) = varrho_z;
Psi(ARz,eps_u_z) = 1;

%------------------------------------
% exogenous process for MP shock(ARM)
%------------------------------------
G0(ARM,u_M) = 1;
G1(ARM,u_M) = varrho_M;
Psi(ARM,eps_u_M) = 1;

%----------------------------
% growth expectation(EXgamma)
%----------------------------
G0(EXgamma,u_gamma) = 1;
G1(EXgamma,xu_gamma) = 1;
Pi(EXgamma,eta_xu_gamma) = 1;

%--------------------------------------------------
% marginal utility of consumption expectation(EXUc)
%--------------------------------------------------
G0(EXUc,Uc) = 1;
G1(EXUc,xUc) = 1;
Pi(EXUc,eta_xUc) = 1;

%----------------------------
% investment expectation(EXi)
%----------------------------
G0(EXi,ivst) = 1;
G1(EXi,xiv) = 1;
Pi(EXi,eta_xiv) = 1;

%------------------------------
% rental rate expectation(EXrk)
%------------------------------
G0(EXrk,rk) = 1;
G1(EXrk,xrk) = 1;
Pi(EXrk,eta_xrk) = 1;

%---------------------------
% Tobin's q expectation(EXq)
%---------------------------
G0(EXq,q) = 1;
G1(EXq,xq) = 1;
Pi(EXq,eta_xq) = 1;

%------------------------------------------
% lagrange multiplier expectation(EXlambda)
%------------------------------------------
G0(EXlambda,lambda) = 1;
G1(EXlambda,xlambda) = 1;
Pi(EXlambda,eta_xlambda) = 1;

%----------------------------------------
% long term bond return expectation(EXRL)
%----------------------------------------
G0(EXRL,RL) = 1;
G1(EXRL,xRL) = 1;
Pi(EXRL, eta_xRL) = 1;

%------------------------------
% inflation expectation(EXinfl)
%------------------------------
G0(EXinfl,infl) = 1;
G1(EXinfl,xinfl) = 1;
Pi(EXinfl,eta_xinfl) = 1;

%----------------------
% wage expectation(EXw)
%----------------------
G0(EXw,w) = 1;
G1(EXw,xw) = 1;
Pi(EXw,eta_xw) = 1;

%-----------------------------------
% lag terms for measurement equation
%-----------------------------------
G0(Lagy,Ly) = 1;
G1(Lagy,y) = 1;

G0(Lagivst,Livst) = 1;
G1(Lagivst,ivst) = 1;

G0(Lagw,Lw) = 1;
G1(Lagw,w) = 1;

G0(Lagb,Lb) = 1;
G1(Lagb,b) = 1;

%-------------------------------
% compute the variance of shocks
%-------------------------------
S = diag([sigma_gamma sigma_U sigma_i sigma_w sigma_p sigma_g sigma_z sigma_M]/100).^2;
%}
end