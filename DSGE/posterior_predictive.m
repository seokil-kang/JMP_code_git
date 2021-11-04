% posterior predictive analysis for monetary policy shock on govt debt valuation

% house keeping
clear;close all;clc;

% choose the version of model
cd benchmark

% length of horizon for impulse response function
H = 100;

% confidence level
conf_lev = .9;

% load posterior sampler
load result/SMC_posterior

% the final stage particle is the posterior sampler
theta = theta(:,:,end);

% add path for submodules
addpath('../functions','../gensys','model','../prior')

if exist('result/posterior_predictive.mat')
    load result/posterior_predictive.mat
else    
    % compute IRF and Debt valuation decomposition
    [IRF DVD FEVD] = IRF_and_DVD(theta,H);
    save('result/posterior_predictive.mat','IRF','DVD','FEVD')
end

% call variable indice
variable_listing

% select variables of interest
var_select = [y infl Q mktb lambda R RL s];
IRF = IRF(var_select,:,:);
FEVD = FEVD(var_select,:,:,:);

% show main result
report_results(IRF,DVD,FEVD,var_select,conf_lev)
cd ..