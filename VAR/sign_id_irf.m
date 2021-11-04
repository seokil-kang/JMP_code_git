% Estimate IRF with sign restriction

% house keeping
clear;clc;close all;

% addpath for required functions
addpath('functions')

% Length of IRF horizon
H = 40;

% load the reduced form estimation result
load HBVAR/result/BVAR_posterior

% call variable indexing
variable_indexing(BVAR_est.data_index);
if exist('y')
    output = y;
elseif exist('x')
    output = x;
end

% denote the VAR equation of sign restriction
shock_index = R;

% measure the VAR dimension
[M,~,~] = size(BVAR_est.Sigma_sampler);

% length of IRF restricted
H_sign = 4;

% sign restriction basket
sign_irf = zeros(M,H_sign);

% sign restriction: contractionary MP shock

% interest rate rises
sign_irf(R,:) = 1;

% output falls
sign_irf(output,1) = -1;

% price falls
sign_irf(Pi,1) = -1;

% primary surplus falls at the impact
sign_irf(ps,1) = -1;

% market value of debt: indetermined
sign_irf(b,1) = -1;

% bond price falls
sign_irf(RB,1) = -1;

% investment falls
sign_irf(i,1) = -1;

% wage and labor: indetermined
sign_irf(w,1) = -1;
sign_irf(n,1) = -1;

% load the result of counterpart DSGE
cd ../DSGE/benchmark
load result/posterior_predictive_25bps.mat
cd model
variable_listing
var_select = [y ivst w n mktb s infl R RL];
IRF_DSGE = IRF(var_select,:,:);
IRF_DSGE = mean(IRF_DSGE,3);
IRF_DSGE = IRF_DSGE(:,1:min(length(IRF_DSGE),H+1));
DVD_DSGE = DVD;
cd ../../../VAR

% run the sign restriction IRF function
tic;
[IRF DVD FEVD] = sign_restriction(BVAR_est,H,sign_irf,shock_index);
toc;

%%
% report results
conf_lev = .9;
% significance level
sig_lev = (1 - conf_lev) / 2;
% horizon to show
H_show = H;
hrz = linspace(0,H,H+1)';
% line width for plot
lw = 1.5;
lw2 = 2.5;
ftsz1 = 18;
ftsz2 = 16;

y_nm_full = ["output";"consumption";"investment";"wage";"labor";"market value of debt";"primary surplus";"govt spending";"transfer";"inflation";"interest rate";"nominal bond return";"output gap";"primary surplus"];
data_index = [1 3 4 5 6 7 10 11 12];
y_nm = y_nm_full(data_index);


N_quantile = 2;
posterior_band_q = linspace(sig_lev/2,.5,N_quantile);
%posterior_band_q = [sig_lev/2 1/6 .5];
%N_quantile = 3;
%mycol = [.0 .0 .5];
mycol = [0, 83, 154]/255/1.5;
mycol2 = [230, 115, 0]/255;
figure('WindowState','maximized','name','IRF to MP Shock','color','w')
tiledlayout(3,3, 'Padding', 'none', 'TileSpacing', 'compact');
for j = 1:M
    IRFj = squeeze(IRF(j,:,:));
    nexttile
    hold on
    yline(0,'color','#a2191f')    
    for jj = 1:N_quantile
        low_quantile = quantile(IRFj,1-posterior_band_q(jj),2);
        high_quantile = quantile(IRFj,posterior_band_q(jj),2);
        if jj == 1
            fill([hrz', flip(hrz')], [low_quantile' flip(high_quantile')],mycol,'facealpha',1/4,'linestyle','none');
        else
            fill([hrz', flip(hrz')], [low_quantile' flip(high_quantile')],mycol,'facealpha',jj/(N_quantile^1.2),'linestyle','none');
        end
    end    
    p1 = plot(hrz,mean(IRFj,2),'o-','color',mycol,'linewidth',lw,'markersize',5);        
    hold off
    set(gca, 'YGrid', 'on', 'XGrid', 'off')
    ticklabel = get(gca,'TickLabel');
    set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz2,'FontWeight','bold');
    xticks([0:8:min(H,H_show)])
    xlim([0 min(H,H_show)])
    title(y_nm(j),'FontSize',ftsz1,'FontName','Consolas','FontWeight','bold');    
end

%report_results_comparison(IRF,IRF_DSGE,BVAR_est.data_index,conf_lev)