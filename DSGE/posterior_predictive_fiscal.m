% posterior predictive analysis for monetary policy shock on govt debt valuation

% house keeping
clear;close all;clc;

% choose the version of model
cd benchmark

% confidence level
conf_lev = .9;

% load posterior sampler
load result/SMC_posterior

% the final stage particle is the posterior sampler
theta = theta(:,:,end);

% add path for submodules
addpath('../functions','../gensys','model','../prior')

% load the posterior predictives
load result/posterior_predictive.mat

% measure the horizon
H = size(IRF,2)-1;

% call variable indice
variable_listing

% select variables of interest
var_select = [s tau g z];
IRF = IRF(var_select,:,:);
IRFs = squeeze(IRF(1,:,:));
IRFt = squeeze(IRF(2,:,:));
IRFg = squeeze(IRF(3,:,:));
IRFz = squeeze(IRF(4,:,:));
FEVD = FEVD(var_select,:,:,:);

cd ..
%%
% significance level
sig_lev = (1 - conf_lev) / 2;

% horizon for x-axis
hrz = linspace(0,H,H+1)';

% line width for plot
lw = 4;
ftsz = 28;
N_quantile = 1;
posterior_band_q = linspace(sig_lev/2,.5,N_quantile);
mycol = {'#00539a' '#ff832b' '#198038' '#a56eff'};
mycols = [0 83 154; 255 131 43; 25 128 56; 165 110 255]/255;
linsty = {'-' ':' '--' '-.'};
mksty = {'none' 'none' 'none' 'none'};
%mksty = {'x' 'o' 'd' 's'};
figure('WindowState','maximized','name','Fiscal Responses to MP Shock','color','w')
tiledlayout(1,1, 'Padding', 'none', 'TileSpacing', 'compact');
nexttile
yline(0,'color','#a2191f') 
for j = 2:4
    IRFj = squeeze(IRF(j,:,:));
    hold on    
    for jj = 1:N_quantile
        low_quantile = quantile(IRFj,1-posterior_band_q(jj),2);
        high_quantile = quantile(IRFj,posterior_band_q(jj),2);
        if jj == 1
            fill([hrz', flip(hrz')], [low_quantile' flip(high_quantile')],mycols(j-1,:),'facealpha',1/10,'linestyle','none');
        else
            fill([hrz', flip(hrz')], [low_quantile' flip(high_quantile')],mycols(j-1,:),'facealpha',jj/(N_quantile^1.2),'linestyle','none');
        end
    end    
    plot(hrz,mean(IRFj,2),'color',mycol{j-1},'linestyle',linsty{j-1},'marker',mksty{j-1},'linewidth',lw);
    hold off
    set(gca, 'YGrid', 'on', 'XGrid', 'off')    
    ticklabel = get(gca,'TickLabel');
    set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz,'FontWeight','bold');
    xlim([0 H])
end
%xlabel('quarter','FontSize',ftsz,'FontName','Consolas','FontWeight','bold');
%ylabel('%','FontSize',ftsz,'FontName','Consolas','FontWeight','bold');
if N_quantile == 1
    legend('','','Tax revenue','','Spending','','Transfer','location','northwest','orientation','vertical','FontSize',ftsz,'FontName','Consolas','FontWeight','bold');
else
    legend('','','','Tax revenue','','','Spending','','','Transfer','location','northwest','orientation','vertical','FontSize',ftsz,'FontName','Consolas','FontWeight','bold');
end
legend boxoff