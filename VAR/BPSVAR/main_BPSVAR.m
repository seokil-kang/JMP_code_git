%-----------------
% HOUSE KEEPING
%----------------

clear all;  
clc; 
close all;

addpath('./auxfiles')
addpath('./results')

% baskets
A_sampler = [];
IMPACT = [];
A0s = [];
Sigma_sampler = [];

%------------------------------------------------------------
% SETTINGS
%------------------------------------------------------------
p = 4;                                  % Number of lags
nex_ = 1;                               % Constant
Horizon = 1;                            % Horizon for calculation of impulse responses
MP = 0;                                 % Tightness of the Minnesota Prior (0: loose; 1: tight)
T0 = p;                                 % Length of pre-sample for Minnesota Prior
nd = 2e4;                               % Number of draws in MC chain
bburn = 0.1*nd;                         % Burn-in period
fflagFEVD = 0;                          % Compute FEVD
differenced = 1;
N_sys = min(nd-bburn,2e4);
N_sys_index = round(linspace(1,nd-bburn,N_sys));

ptileVEC = [0.05 0.16 0.50 0.86 0.95]; % Percentiles

%------------------------------------------------------------
% PRIOR SELECTION -- SIGMA_NU
% IG - Inverse Gamma
% FI - Truncated at pr_trunc x std(M_t)
%------------------------------------------------------------
prior_type = 'FI';
pr_trunc = .5;                      % only used for FIXED

mu0 = 1.00;                             % prior mean
V0 = 0.5^2;                             % prior variance

s0 = 0.02;
nu0 = 2.00;
%------------------------------------------------------------
% Sampler Settings
%------------------------------------------------------------
% this is the mixture probability for the "RW" IG propsal on SIGMA
rwmh_sigma_prob = 0.0;
rwmh_df = 5;                            % Tune-up parameter for mixture proposal distribution for (\Phi,\Sigma)
nu = rwmh_df;                           % Tune-up parameter for mixture proposal distribution for (\Phi,\Sigma)

%---------------------------------------------------------------------
% LOAD DATA
%---------------------------------------------------------------------

i_var_transf =  {};     % Variable that require additional transformations
nCalc = length(i_var_transf);
nlags_ = p;

%-------------------------------------------
% Load data and construst dummy observations
%-------------------------------------------
% load observation
cd ../../data
[YY dateq data_index diff_index] = data_importing([1 3 4 5 6 7 10 11 12],'1960-1','2006-4',differenced);
cd ../VAR/BPSVAR

% load proxy
% load the instrumental variables from DSGE estimates
load ../../DSGE/benchmark/result/estimated_MP_shock
mm = MP_shock(T0+1:end);

nv      = size(YY,2);     %* number of variables */
nobs    = size(YY,1)-T0;  %* number of observations */

if MP == 0
    tau	   =   0.5; 
    d	   =   1;
    w	   =   1;
    lambda =   0.5;	 
    mu	   =   0.5;	 
elseif MP == 1
    tau	   =   0.2; 
    d	   =   3.5;
    w	   =   1;
    lambda =   0.1;	 
    mu	   =   0.1;	 
end

YY0     =   YY(1:T0,:);  
ybar    =   mean(YY0)';
sbar    =   std(YY0)' ;
premom  =   [ybar sbar];

% Generate matrices with dummy observations
hyp = [tau; d; w; lambda; mu];
[YYdum, XXdum, breakss] = varprior_h(nv,nlags_,nex_,hyp,premom);

% Actual observations

YYact = YY(T0+1:T0+nobs,:);
XXact = zeros(nobs,nv*nlags_);

i = 1;

while (i <= nlags_)
    XXact(:,(i-1)*nv+1:i*nv) = YY(T0-(i-1):T0+nobs-i,:);
    i = i+1;
end

XXact = [XXact ones(nobs,1)];

Mobs       = size(mm,1);            
MM         = [ones(Mobs,1) mm]; % Add constant to matrix of proxies

% For uninformative prior uncomment the following two lines
% XXdum = [];
% YYdum = [];

%-------------------------------------------
% Declare objects for estimation
%-------------------------------------------

n = nv;
nIV = size(mm,2);
nshocks=nIV;

e = eye(n); % create identity matrix
aalpha_index = 2:n;
ddelta_index = 1;
a = cell(p,1);

% Define matrices to compute IRFs      
J = [eye(n);repmat(zeros(n),p-1,1)]; % Page 12 RWZ
F = zeros(n*p,n*p);    % Matrix for Companion Form
I  = eye(n);
for i=1:p-1
    F(i*n+1:(i+1)*n,(i-1)*n+1:i*n) = I;
end

% Estimation Preliminaries
X = [XXdum; XXact];
Y = [YYdum; YYact];
T = size(X, 1);
ndum = size(XXdum, 1);
nex = nex_;

% Compute OLS estimates
B = (X'*X)\(X'*Y); % Point estimates
U = Y-X*B;         % Residuals
Sigmau = U'*U/(T-p*n-1);   % Covariance matrix of residuals
F(1:n,1:n*p)    = B(1:n*p,:)';

%------------------------------------------------------------
% MCMC Algorithm
%------------------------------------------------------------

% set preliminaries for priors
N0=zeros(size(X',1),size(X,2));
nnu0=0;
nnuT = T - n*p -nex_;%# +nnu0;
NT = N0 + X'*X;    
Bbar0=B;
S0=Sigmau;
BbarT = NT\(N0*Bbar0 + (X'*X)*B);
ST = (nnu0/nnuT)*S0 + (T/nnuT)*Sigmau + (1/nnuT)*((B-Bbar0)')*N0*(NT\eye(n*p+nex))*(X'*X)*(B-Bbar0); %% Constant (check)
STinv = ST\eye(n);

record=0;     
counter = 0;
fflagEXP = 0;

disp('                                                                  ');
disp('        BAYESIAN ESTIMATION OF VAR VIA BLOCK MCMC                 ');
disp('                                                                  ');

% Drop constant from M
MM = MM(:, 2:end);

bet = 0.01;
signu = 0.04;

%[Q, ~] = qr(randn(n, nIV));
Xstar = randn(n, 1);
Q = Xstar / norm(Xstar);

R=mvnrnd(zeros(n,1),STinv/nnuT,nnuT)';
Sigmadraw=(R*R')\eye(n);
bbeta = B(:);
SigmaB = kron(Sigmadraw,NT\eye(n*p+nex_));
SigmaB = (SigmaB+SigmaB')/2;
Bdraw = mvnrnd(bbeta,SigmaB);

Bdraw= reshape(Bdraw,n*p+nex,n); % Reshape Bdraw from vector to matrix
Udraw = Y-X*Bdraw;      % Store residuals for IV regressions
LC =chol(Sigmadraw,'lower');
A0chol = (LC')\eye(size(LC,1));

lnp0 = loglik_m_given_y(MM, Udraw(ndum+1:end, :), LC, Q, bet, signu);

% Initialize Omega1 for IRFs
Omega1 = [LC(:,nshocks);zeros((p-1)*n,nshocks)];

nq = 1;
acpt_rf   = 0;
acpt_Q    = 0;
Fstar = F;
Sigmadrawstar = Sigmadraw;

% MCMC Chain 
% Define objects that store the draws
Ltilde = zeros(nd-bburn,Horizon+1,n,nshocks);                      % define array to store IRF
LtildeAdd = zeros(nd-bburn,Horizon+1,n+nCalc,nshocks);
irfCalc = zeros(nd-bburn,Horizon+1,nCalc,nshocks);                     % store labor productivity IRF
W = zeros(nd-bburn,Horizon+1,n+nCalc,nshocks);                         % define array to store FVD
EETA  = zeros(nd-bburn,n-1);
ppsi_levels = zeros(nd-bburn,n);
ppsi_diff   = zeros(nd-bburn,n-1);
REL = zeros(nd-bburn,1);
BET = zeros(nd-bburn,1);
SIG = zeros(nd-bburn,1);
A0MAT = zeros(nd-bburn, n, n);
ApMAT = zeros(nd-bburn, n*p+nex, n);
while record<nd
    
    %------------------------------------------------------------
    % Gibbs Sampling Algorithm
    %------------------------------------------------------------
    % STEP ONE: Draw from the B, SigmaB | Y
    %------------------------------------------------------------
    % Step 1: Draw from the marginal posterior for Sigmau p(Sigmau|Y,X)
    if (rand()<rwmh_sigma_prob)     % draw sigma* ~ IW(sigma, nu), phi*|sigma*, Y
        nu = nobs;
        Sigmadrawstar = iwishrnd(nu*Sigmadraw, nu+n+1);
        
        bbeta = B(:);
        SigmaB = kron(Sigmadrawstar,NT\eye(n*p+nex));
        SigmaB = (SigmaB+SigmaB')/2;
        Bdrawstar = mvnrnd(bbeta,SigmaB);
        
        Bdrawstar= reshape(Bdrawstar,n*p+nex,n);% Reshape Bdraw from vector to matrix
        Udrawstar = Y-X*Bdrawstar;              % Store residuals for IV regressions
        
        LCstar    = chol(Sigmadrawstar,'lower');
        A0cholstar    = (LCstar')\eye(size(LCstar,1));
        
        Fstar(1:n,1:n*p)    = Bdrawstar(1:n*p,:)';
        
        lnp1 = loglik_m_given_y(MM, Udrawstar(ndum+1:end, :), LCstar, Q(:,nIV), bet, signu);
        
        q1 = pdf_ln_iwish(nu*Sigmadraw, nu+n+1, Sigmadrawstar);
        q0 = pdf_ln_iwish(nu*Sigmadrawstar, nu+n+1, Sigmadraw);
        
        psigma1 = pdf_ln_iwish(ST, nnuT, Sigmadrawstar);
        psigma0 = pdf_ln_iwish(ST, nnuT, Sigmadraw);            
        
        % metropolis hastings
        alp = exp((lnp1+psigma1) - (lnp0+psigma0) - (q1 - q0));
        
        
    else                            % draw phi*, sigma* ~ phi, sigma | Y
        
        R=mvnrnd(zeros(n,1),STinv/nnuT,nnuT)';
        Sigmadrawstar=(R*R')\eye(n);
        
        % Step 2: Taking newSigma as given draw for B using a multivariate normal    
        bbeta = B(:);
        SigmaB = kron(Sigmadrawstar,NT\eye(n*p+nex));
        SigmaB = (SigmaB+SigmaB')/2;
        Bdrawstar = mvnrnd(bbeta,SigmaB);
        % Storing unrestricted draws
        
        Bdrawstar= reshape(Bdrawstar,n*p+nex,n);% Reshape Bdraw from vector to matrix
        Udrawstar = Y-X*Bdrawstar;      % Store residuals for IV regressions
        
        LCstar    = chol(Sigmadrawstar,'lower');
        A0cholstar    = (LCstar')\eye(size(LCstar,1));
        
        Fstar(1:n,1:n*p)    = Bdrawstar(1:n*p,:)';
        
        lnp1 = loglik_m_given_y(MM, Udrawstar(ndum+1:end, :), LCstar, Q(:,nIV), bet, signu);
        
        % metropolis hastings
        alp = exp(lnp1 - lnp0);
    end
    
    if rand < alp
        lnp0 = lnp1;
        
        Udraw = Udrawstar;
        LC = LCstar;
        A0chol = A0cholstar;
        
        F = Fstar;
        Sigmadraw = Sigmadrawstar;
        Bdraw = Bdrawstar;           
        acpt_rf = acpt_rf+1;
    end
    
    %------------------------------------------------------------
    % STEP TWO: Draw from Q distribution
    %------------------------------------------------------------
    %            [Qstar, ~] = qr(randn(n, nIV));
    Xstar = randn(n, 1);
    Qstar = Xstar / norm(Xstar);
    
    lnp1 = loglik_m_given_y(MM, Udraw(ndum+1:end, :), LC, Qstar(:,nIV), bet, signu);
    
    % metropolis hastings
    alp = exp(lnp1 - lnp0);
    if rand < alp
        lnp0 = lnp1;
        Q = Qstar;
        acpt_Q = acpt_Q+1;
    end
    
    % normalize A0
    A0 = A0chol*Q;
    if A0(1) < 0
        Q = -Q;
    end    
    
    %------------------------------------------------------------
    % STEP THREE: BET and SIGMA from MVN-IG / MVN-IW
    %------------------------------------------------------------
    scale_mat = Q' / LC;
    Xe = (scale_mat*Udraw(ndum+1:end, :)')';
    
    % bayesian linear regresion model
    Bhat = inv(Xe'*Xe)*Xe'*MM;
    Vp = inv(Xe'*Xe + inv(V0/signu^2));
    mup = Vp * (Xe'*Xe*Bhat + inv(V0)*mu0);
    
    if (prior_type == 'IG')          
        nu1 = size(Xe, 1)-1 + nu0;
        s1 = ( (MM-Xe*mup)'*(MM-Xe*mup) + nu0*s0^2)/nu1;
        
        signu = igrand(s1, nu1);
    elseif (prior_type == 'FI')     
        signu = pr_trunc*std(MM);
    end
    
    bet = mvnrnd(mup, signu^2*Vp);
    
    lnp0 = loglik_m_given_y(MM, Udrawstar(ndum+1:end, :), LCstar, Q(:,nIV), bet, signu);
    
    record=record+1;
    counter = counter +1;
    if counter==0.05*nd
        disp(['         DRAW NUMBER:   ', num2str(record)]);
        disp('                                                                  ');
        disp(['     REMAINING DRAWS:   ', num2str(nd-record)]);
        disp('                                                                  ');
        counter = 0;
        
        fprintf('RF acpt rate: %5.3f\n',   acpt_rf/record)
        fprintf('Q acpt rate: %5.3f\n',    acpt_Q/(record*nq))
    end
    
    if record > bburn
        
        ffactor = LC*Q;
        A0 = A0chol*Q;
        Aplus = Bdraw(1:n*p,:)*A0;
        
        % Compute cumulative coefficients 
        a0 = A0(:,1);
        for l=1:p
            a{l} = Aplus((l-1)*n+1:l*n,1);
        end
        aalpha = zeros(p+1,n-1);
        for j=1:n-1
            jj=aalpha_index(1,j);
            aalpha(1,j)  = -((e(:,jj)'*a0));
            for l=1:p
                aalpha(l+1,j)  = ((e(:,jj)'*a{l}));
            end
        end
        ddelta = zeros(p+1,1);
        ddelta(1,1) = ((e(:,ddelta_index )'*a0));
        for l=1:p
            ddelta(l+1,1)  = -((e(:,ddelta_index )'*a{l}));
        end
        tmp_levels   = zeros(1,n-1);
        
        for l=0:p
            tmp_levels = aalpha(l+1,1:n-1) + tmp_levels;
        end
        % variables in levels
        tmp_diff   = zeros(1,n-1);
        for l=0:p
            for ii=0:l
                tmp_diff = aalpha(ii+1,:) +  tmp_diff ;
            end
        end
        tmpR=0;
        for l=0:p
            tmpR=ddelta(l+1,1)+tmpR;
        end
        
        aalphaEta = aalpha./a0(1);
        ddeltaEta = ddelta./a0(1);
        den = ddeltaEta(1) + sum(ddeltaEta(2:end));       
        num = sum(aalphaEta);
        numtemp = sum(cumsum(aalphaEta));
        REL(record-bburn,1) = bet^2/(bet^2 + signu^2);
        EETA(record-bburn,:) = -A0(2:n,1)./A0(1,1);
        ppsi_levels(record-bburn,:) = [sum(ddeltaEta(2:end)) num];
        ppsi_diff(record-bburn,1:n-1)   = numtemp;
        BET(record-bburn,1) = bet;
        SIG(record-bburn,1) = signu;
        IRF_T    = vm_irf(F,J,ffactor,Horizon+1,n,Omega1);        
        Ltilde(record-bburn,:,:,:) = IRF_T(1:Horizon+1,:,:);
        
        if nCalc
            for ii = 1:length(i_transf)
                irfCalc(record-bburn,:,ii,:) = cumsum(squeeze(IRF_T(:,i_transf(ii),:)));
            end
        end
        
        if fflagFEVD ==1
            W(record-bburn,:,:,:)=variancedecompositionFD(F,J,Sigmadraw,ffactor,n,Horizon,i_transf);
        end
        %%%%%%%
        A_sampler = cat(3,A_sampler,Bdraw);
        Sigma_sampler = cat(3,Sigma_sampler,Sigmadraw);
        IMPACT = [IMPACT ffactor];
        A0s = [A0s A0];
        %%%%%%%
    end
end

LtildeAdd(:,:,1:n,:) = Ltilde;
LtildeAdd(:,:,n+1:n+nCalc,:) = irfCalc;
SVAR.LtildeImpact = LtildeAdd(:,1,:,:);
LtildeFull = quantile(LtildeAdd,ptileVEC);
SVAR.LtildeFull = permute(LtildeFull,[3,2,1,4]);
WhFull = quantile(W,ptileVEC);
SVAR.WhFull = permute(WhFull,[3,2,1,4]);
SVAR.EETAFull = quantile(EETA,ptileVEC);
SVAR.PPSIFull = quantile(ppsi_levels,ptileVEC);
SVAR.PPSIDFull = quantile(ppsi_diff,ptileVEC);
SVAR.RELFull  = quantile(REL,ptileVEC);
SVAR.BETFull  = quantile(BET,ptileVEC);
SVAR.SIGFull  = quantile(SIG,ptileVEC);
SVAR.SIG = SIG;
SVAR.BET = BET;
SVAR.acpt_rf = acpt_rf/record;
SVAR.acpt_Q  = acpt_Q/(record*nq);
SVAR.p = p;
SVAR.nd = nd;
SVAR.bburn = bburn;
SVAR.fflagFEVD = fflagFEVD;
SVAR.pr_trunc = 0;
SVAR.pr_truncFlag = 0;
SVAR.EETA = EETA;
%%
BVAR_est.A_sampler = A_sampler(:,:,N_sys_index);
BVAR_est.Sigma_sampler = Sigma_sampler(:,:,N_sys_index);
BVAR_est.impact = IMPACT(:,N_sys_index);
BVAR_est.data_index = data_index;
BVAR_est.diff_index = diff_index;

% save the posterior sampler
save('results/BVAR_posterior.mat','BVAR_est')