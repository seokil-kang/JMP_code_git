% extract the estimated monetary policy shock from the DSGE

% house keeping
clear;close all;clc;

% choose the version of model
cd benchmark

% confidence level
conf_lev = .9;

% load posterior sampler
load result/SMC_posterior

% posterior sampler is the particles of the last stage
theta = theta(:,:,end);

% measure sizes
[n_param N] = size(theta);

% add path for submodules
addpath('../../data','../functions','../gensys','model','../prior','../sims_minwel')

% prior info setup
[prior_info nm]= prior_info_setup;

% call variable indice
variable_listing

% measurement error scale
ME_scale = 1;

% MP shock basket
eps_M = [];

% extracting
for j = 1:N
    
    % pick a parameter set
    thetaj = theta(:,j);
    
    % kalman smoothing
    [X_smooth X_update G C M] = kalman_smoother(thetaj,y_data,ME_scale);
    
    % AR(1) MP shock process
    u_Mj = X_smooth(u_M,:)';
    
    % AR(1) MP shock parameters
    varrho_M = thetaj(find(nm == 'varrho_M'));
    sigma_M = thetaj(find(nm == 'sigma_M'))/1e2;
    
    % structural shock from standard normal
    eps_M = [eps_M (u_Mj(2:end) - varrho_M*u_Mj(1:end-1)) / sigma_M];
end
%% reporting plot
% significance level
sig_lev = (1 - conf_lev) / 2;

% credible set of the shock
eps_MQ = quantile(eps_M,[.5 sig_lev 1-sig_lev],2);

% load Romer and Romer shock
RR = readmatrix('RomerandRomerDataAppendix.xls','sheet','DATA BY MONTH','range','B2:B373');
dateRR = [datetime(1966,1,1,'format','yyyy-Q'):calmonths(1):datetime(1996,12,31,'format','yyyy-Q')]';

% quarterize it
RRQ = sum(reshape(RR,3,length(RR)/3))';
shock_corr = mean(corr(RRQ,eps_M(find(dateq==dateRR(1)):find(dateq <= dateRR(end),1,'last'),:)));
% plot the shock
N_quantile = 2;
posterior_band_q = linspace(sig_lev/2,.5,N_quantile);
mycol = {'#00539a' '#ff832b' '#198038' '#a56eff'};
mycols = [0 83 154; 255 131 43; 25 128 56; 165 110 255]/255;
linsty = {'-' '-.' ':' '--'};
mksty = {'none' 'none' 'none' 'none'};
%mksty = {'o' 's' 'd' 'x'};
lw = 3;
   
figure('name','Estimated Monetary Policy Shock','WindowState','maximized','color','w')
tiledlayout(1,1, 'Padding', 'none', 'TileSpacing', 'compact');
nexttile
hold on
fill([dateq(3:end)', flip(dateq(3:end)')], [eps_MQ(3:end,2)' flip(eps_MQ(3:end,3)')],mycols(1,:),'facealpha',1/4,'linestyle','none');
p1 = plot(dateq(3:end),eps_MQ(3:end,1),'linewidth',lw,'color',mycol{1},'linestyle',linsty{1},'marker',mksty{1});
p2 = plot(dateRR(1:3:end),RRQ,'linewidth',lw,'color',mycol{2},'linestyle',linsty{2},'marker',mksty{2});
hold off
set(gca,'YGrid','on')
recessionplot
set(gca, 'YGrid', 'on', 'XGrid', 'off')
ticklabel = get(gca,'TickLabel');
set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',24,'FontWeight','bold');
legend([p1 p2],'Estimated MP shock', 'Romer&Romer', 'location','northwest','orientation','vertical','FontSize',24,'FontName','Consolas','FontWeight','bold');
legend boxoff
title(sprintf('correlation between two shocks = %.2f',shock_corr),'FontName','Consolas','fontsize',24,'FontWeight','bold');

%% save result
dateshock = dateq;
MP_shock = eps_MQ(:,1);

%save('result/estimated_MP_shock','dateshock','MP_shock')

cd ..