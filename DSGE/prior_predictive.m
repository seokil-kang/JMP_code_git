% prior predictive analysis for monetary policy shock on govt debt valuation

% house keeping
clear;close all;clc;

% choose the version of model
cd benchmark

% # of prior draws
N = 2e4;

% length of horizon for impulse response function
H = 40;

% confidence level
conf_lev = .9;

% add path for submodules
addpath('../functions','../gensys','model','../prior')

if exist('result/prior_predictive.mat')
    load result/prior_predictive.mat
else    
    % prior info setup
    [prior_info nm]= prior_info_setup;
    
    % draw from parameter distribution
    theta = prior_drawing(prior_info,N);
    
    % compute IRF and Debt valuation decomposition
    [IRF DVD FEVD] = IRF_and_DVD(theta,H);
    save('result/prior_predictive.mat','IRF','DVD','FEVD','-v7.3')
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