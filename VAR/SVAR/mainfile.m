% SVAR estimation: Arias, Caldara, Rubio-Ramirez, JME, 2019

%==============
% house keeping
%==============
clear all;close all; clc;

differenced = 0;

%=============
% data loading
%=============
cd ../../data
[y_obs dateq data_index diff_index] = data_importing([1 3 4 5 6 7 10 11 12],'1960-1','2006-4',differenced);
cd ../VAR/SVAR

%=================
% estimation setup
%=================
rng('default'); % reinitialize the random number generator to its startup configuration
rng(0);         % set seed
addpath('helpfunctions'); % set path to helper functions
nlag = 4; % number of lags
nvar = size(y_obs,2); % number of endogenous variables
nex = 1; % set equal to 1 if a constant is included; 0 otherwise
m = nvar*nlag + nex; % number of exogenous variables
nd        = 2e4; % number of orthogonal-reduced-form (B,Sigma,Q) draws
iter_show = 1e3; % display iteration every iter_show draws
horizon   = 0; % maximum horizon for IRFs
%index     = 40; % define  horizons for the FEVD
horizons  = 0:1; % horizons to restrict
NS        = 1 + numel(horizons); % number of objects in F(A_{0},A_{+}) to which we impose sign and zero restrictios: F(THETA)=[A_{0};L_{0},...,L_{horizons}]
e         = eye(nvar); % create identity matrix
maxdraws  = 2e4; % max number of importance sampling draws
conjugate = 'structural'; % structural or irfs or empty

%==========================================================================
%% identification: declare Ss and Zs matrices
%==========================================================================
% restrictions on A0 and/or IRFs
variable_indexing(data_index)
if exist('y')
    output = y;
else
    output = x;
end

% sign restrictions
S = cell(nvar,1);
for ii=1:nvar
    S{ii}=zeros(0,nvar*NS);
end
ns1  = 5;
S{1} = zeros(ns1,nvar*NS);
S{1}(1,output)   = -1;
S{1}(2,Pi)   = -1;
S{1}(3,R)   =  1;
S{1}(4,nvar+R)   =  1;
S{1}(5,nvar+RB)   =  -1;
%S{1}(6,nvar+ps)   =  -1;

% zero restrictions
Z=cell(nvar,1);
for i=1:nvar
    Z{i}=zeros(0,nvar*NS);
end

nz1  = 6;
Z{1} = zeros(nz1,nvar*NS);
Z{1}(1,i) = 1;
Z{1}(2,w) = 1;
Z{1}(3,n) = 1;
Z{1}(4,b) = 1;
Z{1}(5,ps) = 1;
Z{1}(6,RB) = 1;

%==========================================================================
%% Setup info
%==========================================================================
info=SetupInfo(nvar,m,Z,@(x)chol(x));

% ZF(A_{0},A_{+})
info.nlag     = nlag;
info.horizons = horizons;
info.ZF       = @(x,y)ZF(x,y);

% functions useful to compute the importance sampler weights
iw_info = info;
fs      = @(x)ff_h(x,iw_info);
r       = @(x)ZeroRestrictions(x,iw_info);

if strcmp(conjugate,'irfs')==1
    fo              = @(x)f_h(x,iw_info);
    fo_str2irfs     = @(x)StructuralToIRF(x,iw_info);
    fo_str2irfs_inv = @(x)IRFToStructural(x,iw_info);
    r_irfs          = @(x)IRFRestrictions_more_general(x,iw_info); 
end


% function useful to check the sign restrictions
fh_S_restrictions  = @(x)SF(x,iw_info,S);

%==========================================================================
%% write data in Rubio, Waggoner, and Zha (RES 2010)'s notation
%==========================================================================
% yt(t) A0 = xt(t) Aplus + constant + et(t) for t=1...,T;
% yt(t)    = xt(t) B     + ut(t)            for t=1...,T;
% x(t)     = [yt(t-1), ... , yt(t-nlag), constant];
% matrix notation yt = xt*B + ut;
% xt=[yt_{-1} ones(T,1)];
yt = y_obs(nlag+1:end,:);
T  = size(yt,1);
xt = zeros(T,nvar*nlag+nex);
for i=1:nlag
    xt(:,nvar*(i-1)+1:nvar*i) = y_obs((nlag-(i-1)):end-i,:) ;
end
if nex>=1
    xt(:,nvar*nlag+nex)=ones(T,1);
end
% write data in Zellner (1971, pp 224-227) notation
Y = yt; % T by nvar matrix of observations
X = xt; % T by (nvar*nlag+1) matrix of regressors


%% prior for reduced-form parameters
nnuBar              = 0;
OomegaBarInverse    = zeros(m);
PpsiBar             = zeros(m,nvar);
PpsiBar             = zeros(m,nvar) + [eye(nvar); zeros(m-nvar,nvar)]; % MN RW prior
PphiBar             = zeros(nvar);

%% posterior for reduced-form parameters
nnuTilde            = T +nnuBar;
OomegaTilde         = (X'*X  + OomegaBarInverse)\eye(m);
OomegaTildeInverse  =  X'*X  + OomegaBarInverse;
PpsiTilde           = OomegaTilde*(X'*Y + OomegaBarInverse*PpsiBar);
PphiTilde           = Y'*Y + PphiBar + PpsiBar'*OomegaBarInverse*PpsiBar - PpsiTilde'*OomegaTildeInverse*PpsiTilde;
PphiTilde           = (PphiTilde'+PphiTilde)*0.5;


%% useful definitions
% definitios used to store orthogonal-reduced-form draws, volume elements, and unnormalized weights
Bdraws         = cell([nd,1]); % reduced-form lag parameters
Sigmadraws     = cell([nd,1]); % reduced-form covariance matrices
Qdraws         = cell([nd,1]); % orthogonal matrices
storevefh      = zeros(nd,1);  % volume element f_{h}
storevegfhZ    = zeros(nd,1);  % volume element g o f_{h}|Z
uw             = zeros(nd,1);  % unnormalized importance sampler weights

if strcmp(conjugate,'irfs')==1
    storevephi      = zeros(nd,1);  % volume element f_{h}
    storevegphiZ    = zeros(nd,1);  % volume element g o f_{h}|Z
end

% definitions related to IRFs; based on page 12 of Rubio, Waggoner, and Zha (RES 2010)
J      = [e;repmat(zeros(nvar),nlag-1,1)];
A      = cell(nlag,1);
extraF = repmat(zeros(nvar),1,nlag-1);
F      = zeros(nlag*nvar,nlag*nvar);
for l=1:nlag-1
    F((l-1)*nvar+1:l*nvar,nvar+1:nlag*nvar)=[repmat(zeros(nvar),1,l-1) e repmat(zeros(nvar),1,nlag-(l+1))];
end

% definition to facilitate the draws from B|Sigma
hh              = info.h;
cholOomegaTilde = hh(OomegaTilde)'; % this matrix is used to draw B|Sigma below


%% initialize counters to track the state of the computations
counter = 1;
record  = 1;
count   = 0;
tStart = tic;
while record<=nd
    
    
    %% step 1 in Algorithm 2
    Sigmadraw     = iwishrnd(PphiTilde,nnuTilde);
    cholSigmadraw = hh(Sigmadraw)';
    Bdraw         = kron(cholSigmadraw,cholOomegaTilde)*randn(m*nvar,1) + reshape(PpsiTilde,nvar*m,1);
    Bdraw         = reshape(Bdraw,nvar*nlag+nex,nvar);
    % store reduced-form draws
    Bdraws{record,1}     = Bdraw;
    Sigmadraws{record,1} = Sigmadraw;
    
   
    %% steps 2:4 of Algorithm 2
    w           = DrawW(iw_info);   
    x           = [vec(Bdraw); vec(Sigmadraw); w];
    structpara  = ff_h_inv(x,iw_info);
    
    % store the matrix Q associated with step 3
    Qdraw            = SpheresToQ(w,iw_info,Bdraw,Sigmadraw);
    Qdraws{record,1} = reshape(Qdraw,nvar,nvar);
    
   
    %% check if sign restrictions hold
    signs      = fh_S_restrictions(structpara);
    
    tmpA0 = reshape(structpara(1:nvar*nvar),nvar,nvar);
    tmpA0inv = inv(tmpA0);
    impact_RB = tmpA0inv(1,9);
    impact_PS = tmpA0inv(1,6);
    ppsiy = -tmpA0(1,1)/tmpA0(8,1);
    ppsip = -tmpA0(7,1)/tmpA0(8,1);
    rhoR = -structpara(nvar*nvar+8)/tmpA0(8,1);
    
    if (sum(signs{1}*e(:,1)>0))==size(signs{1}*e(:,1),1)
    %if (sum(signs{1}*e(:,1)>0))==size(signs{1}*e(:,1),1) && ppsiy>0 && ppsiy<1e5 && ppsip>0 && ppsip<1e5
        
        count=count+1;
        
        %% compute importance sampling weights
       
        switch conjugate
            
            case 'structural'
                
                
                storevefh(record,1)   = (nvar*(nvar+1)/2)*log(2)-(2*nvar+m+1)*LogAbsDet(reshape(structpara(1:nvar*nvar),nvar,nvar));
                storevegfhZ(record,1) = LogVolumeElement(fs,structpara,r); 
                uw(record,1)          = exp(storevefh(record,1) - storevegfhZ(record,1));
                
            case 'irfs'
                
                irfpara                = fo_str2irfs(structpara);
                storevephi(record,1)   = LogVolumeElement(fo,structpara)   + LogVolumeElement(fo_str2irfs_inv,irfpara);
                storevegphiZ(record,1) = LogVolumeElement(fs,structpara,r) + LogVolumeElement(fo_str2irfs_inv,irfpara,r_irfs); 
                uw(record,1)           = exp(storevephi(record,1) - storevegphiZ(record,1));
                
            otherwise
                
                uw(record,1) = 1;
                
        end
        
    else
        
        uw(record,1) = 0;
        
    end
    
    if counter==iter_show
        
        display(['Number of draws = ',num2str(record)])
        display(['Remaining draws = ',num2str(nd-(record))])
        counter =0;
        
    end
    counter = counter + 1;
    record=record+1;
    
end


tElapsed = toc(tStart)
imp_w  = uw/sum(uw);
ne = floor(1/sum(imp_w.^2));


%% store draws
Ltilde        = zeros(horizon+1,nvar,nvar,ne); % define array to store IRF
A0tilde       = zeros(nvar,nvar,ne); % define array to store A0
Aplustilde    = zeros(m,nvar,ne); % define array to store Aplus
hist_is_draws = zeros(ne,1);   % define array to store draws from importance sampler
Sigma_sampler = zeros(nvar,nvar,ne);
Q_sampler = zeros(nvar,nvar,ne);


for s=1:min(ne,maxdraws)
    
    %% draw: B,Sigma,Q
    is_draw = randsample(1:size(imp_w,1),1,true,imp_w);
    hist_is_draws(s,1)=is_draw;
    Sigma_sampler(:,:,s) = Sigmadraws{is_draw};
    Q_sampler(:,:,s) = Qdraws{is_draw};

    Bdraw       = Bdraws{is_draw,1};
    Sigmadraw   = Sigmadraws{is_draw,1};
    Qdraw       = Qdraws{is_draw,1};
    
    
    x=[reshape(Bdraw,m*nvar,1); reshape(Sigmadraw,nvar*nvar,1); Qdraw(:)];
    structpara = f_h_inv(x,info);
    
    
    LIRF =IRF_horizons(structpara, nvar, nlag, m, 0:horizon);
    
    
    for h=0:horizon
        Ltilde(h+1,:,:,s) =  LIRF(1+h*nvar:(h+1)*nvar,:);
    end
    
    
    A0tilde(:,:,s)    = reshape(structpara(1:nvar*nvar),nvar,nvar);
    Aplustilde(:,:,s) = reshape(structpara(nvar*nvar+1:end),m,nvar);
    
    
end

A0tilde    = A0tilde(:,:,1:s);
Aplustilde = Aplustilde(:,:,1:s);
Ltilde     = Ltilde(:,:,:,1:s);

% note that inv(chol(Sigmadraw,'lower')*Qdraw)' = A0tilde
% where A0tilde goes to the SVAR
% y_t' * A0tilde = y_t-1' * Aplustilde + u_t'

cd results
savefile='results.mat';
save(savefile,'Ltilde','A0tilde','Aplustilde','Sigma_sampler','Q_sampler','imp_w','ne','data_index','diff_index');
cd ..