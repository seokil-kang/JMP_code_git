% hierarchical Bayesian VAR estimation by Giannone, Lenza & Primiceri (2015)

% house keeping
clear; close all; clc;

% # of lags
p = 4;

% # of draws
N = 2e4;

% data differencing
differenced = 0;

% add path for submodules
addpath('subroutines')

% load data
cd ../../data
[y dateq data_index diff_index] = data_importing([1 3 4 5 6 7 10 11 12],'1960-1','2006-4',differenced);
cd ../VAR/HBVAR

% run BVAR
rng(0)
res = bvarGLP(y,p,'mcmc',1,'MCMCconst',.5,'Ndraws',N);

% report the acceptance rate of MCMC
fprintf('MCMC acceptance rate = %.2f%s\n',res.mcmc.ACCrate*100,'%')
accrate = (res.mcmc.lambda(2:end)~=res.mcmc.lambda(1:end-1));

% systematic sampling
N = size(res.mcmc.beta,3);
N_sys = min(N,2e4);
N_sys_index = round(linspace(1,N,N_sys));

% take all the results into a structure form
BVAR_est.A_sampler = res.mcmc.beta(:,:,N_sys_index);
BVAR_est.Sigma_sampler = res.mcmc.sigma(:,:,N_sys_index);
BVAR_est.Y = y;
BVAR_est.dateq = dateq;
BVAR_est.data_index = data_index;
BVAR_est.diff_index = diff_index;
%%
% save the posterior sampler
save('result/BVAR_posterior.mat','BVAR_est')

delete *.dat
delete g1.mat g2.mat g3.mat