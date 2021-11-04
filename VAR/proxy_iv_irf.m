% Estimate IRF with proxy instrument variable

% house keeping
clear;clc;close all;

% addpath for required functions
addpath('functions')

% Length of IRF horizon
H = 40;

% confidence level
conf_lev = .9;

% load the reduced form estimation result
load HBVAR/result/BVAR_posterior

% call variable indexing
variable_indexing(BVAR_est.data_index);

% denote the policy index
policy_index = R;

% measure the VAR dimension size
[M,~,~] = size(BVAR_est.Sigma_sampler);

% measure the degree of freedom and # of draws
[K,~,N] = size(BVAR_est.A_sampler);

% is there a constant?
if mod(K,M) == 0
    regressor_constant = 0;
elseif mod(K,M) == 1
    regressor_constant = 1;
end

% lag length is...
p = (K-regressor_constant) / M;

% construct VAR setup
Y_full = BVAR_est.Y;

% construct the dependent variable
Y = Y_full(p+1:end,:);

% cut the date variable as well
dateq = BVAR_est.dateq(p+1:end);

% load the instrumental variables from DSGE estimates
load ../DSGE/benchmark/result/estimated_MP_shock

% cut some suspicious initial points from Kalman smoother
tcut = 1;
Z = MP_shock(tcut:end);
dateshock = dateshock(tcut:end);

% sample period matching between VAR and its proxy
t1 = max(dateq(1),dateshock(1));
t2 = min(dateq(end),dateshock(end));

% construct the regressor matrix
X = [];
for j = 1:p
    X = [X Y_full(p-j+1:end-j,:)];
end
if regressor_constant == 1
    X = [ones(length(Y_full)-p,1) X];
end

% align the sample and IV
Y = Y(find(dateq == t1):find(dateq == t2),:);
X = X(find(dateq == t1):find(dateq == t2),:);
Z = Z(find(dateshock == t1):find(dateshock == t2),:);

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

% run the proxy IV identification IRF function
[IRF DVD FEVD Phi] = proxy_iv_25bps(Y,X,BVAR_est,H,Z,policy_index);

% report results
%report_results(IRF,DVD,FEVD,BVAR_est.data_index,conf_lev)
report_results_comparison(IRF,IRF_DSGE,BVAR_est.data_index,conf_lev)
%% counter-factuals
dvd = rmoutliers(DVD,'percentiles',[(1-conf_lev)/2 1-(1-conf_lev)/2]*1e2,1);
we = dvd(:,2) - dvd(:,1) - dvd(:,3);
dvd_dsge = rmoutliers(DVD_DSGE,'percentiles',[(1-conf_lev)/2 1-(1-conf_lev)/2]*1e2,1);
pi_cf = dvd(:,1) - (we-dvd(:,4));

Cell{1} = {'DSGE','VAR','WE'};
Dell{1} = {'DSGE','VAR','CF'};
ftsz = 20;
ftsz2 = 24;
mksz = 8;
cpsz = 16;
lw = 3;
lw2 = 2;
aphr = .5;
dvd_title  = ["Inflation at impact" "Bond price revaluation" "Discount rate effect" "Fiscal backing" "Wealth effect" "inflation coun-fact"];
x1 = 1;
x2 = 3;
x3 = 2;
xgp = .5;
figure('WindowState','maximized','name','VAR Estimated Fiscal Backing','color','w')
tiledlayout(1,4, 'Padding', 'none', 'TileSpacing', 'compact');    
nexttile        
hold on
yline(0,'--','color','r','linewidth',1)
j = 1;
plot([x1 x2],[mean(dvd_dsge(:,j)) mean(dvd(:,j))],':','color',[1 1 1]*aphr,'linewidth',lw2);
errorbar(x1,mean(dvd_dsge(:,j)),mean(dvd_dsge(:,j))-min(dvd_dsge(:,j)),max(dvd_dsge(:,j))-mean(dvd_dsge(:,j)),'s','color','#00539a','MarkerSize',mksz,'MarkerEdgeColor','#ff832b','MarkerFaceColor','#ff832b','linewidth',lw,'capsize',cpsz)
errorbar(x2,mean(dvd(:,j)),mean(dvd(:,j))-min(dvd(:,j)),max(dvd(:,j))-mean(dvd(:,j)),'s','color','#00539a','MarkerSize',mksz,'MarkerEdgeColor','#198038','MarkerFaceColor','#198038','linewidth',lw,'capsize',cpsz)
hold off
set(gca, 'YGrid', 'on', 'XGrid', 'off')
ticklabel = get(gca,'TickLabel');
set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz2,'FontWeight','bold');
set(gca, 'XTick',[x1 x2], 'XTickLabel',Dell{1})
xlim([x1-xgp x2+xgp])
title(dvd_title(j),'FontName','Consolas','fontsize',ftsz,'FontWeight','bold');
for j = 2:3
    nexttile
    if j == 3
        yline(0,'--','color','r','linewidth',1)
    end
    hold on
    plot([x1 x2],[mean(dvd_dsge(:,j)) mean(dvd(:,j))],':','color',[1 1 1]*aphr,'linewidth',lw2);
    errorbar(x1,mean(dvd_dsge(:,j)),mean(dvd_dsge(:,j))-min(dvd_dsge(:,j)),max(dvd_dsge(:,j))-mean(dvd_dsge(:,j)),'s','color','#00539a','MarkerSize',mksz,'MarkerEdgeColor','#ff832b','MarkerFaceColor','#ff832b','linewidth',lw,'capsize',cpsz)
    errorbar(x2,mean(dvd(:,j)),mean(dvd(:,j))-min(dvd(:,j)),max(dvd(:,j))-mean(dvd(:,j)),'s','color','#00539a','MarkerSize',mksz,'MarkerEdgeColor','#198038','MarkerFaceColor','#198038','linewidth',lw,'capsize',cpsz)
    hold off
    set(gca, 'YGrid', 'on', 'XGrid', 'off')
    ticklabel = get(gca,'TickLabel');
    set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz2,'FontWeight','bold');
    set(gca, 'XTick',[x1 x2], 'XTickLabel',Cell{1})
    xlim([x1-xgp x2+xgp])
    if j == 3
        ylim([-11.5 2.5])
    end
    title(dvd_title(j),'FontName','Consolas','fontsize',ftsz,'FontWeight','bold');
end
nexttile
hold on
yline(0,'--','color','r','linewidth',1)
j = 4;
plot([x1 x3 x2],[mean(dvd_dsge(:,j)) mean(dvd(:,j)) mean(we)],':','color',[1 1 1]*aphr,'linewidth',lw2);
errorbar(x1,mean(dvd_dsge(:,j)),mean(dvd_dsge(:,j))-min(dvd_dsge(:,j)),max(dvd_dsge(:,j))-mean(dvd_dsge(:,j)),'s','color','#00539a','MarkerSize',mksz,'MarkerEdgeColor','#ff832b','MarkerFaceColor','#ff832b','linewidth',lw,'capsize',cpsz)
errorbar(x3,mean(dvd(:,j)),mean(dvd(:,j))-min(dvd(:,j)),max(dvd(:,j))-mean(dvd(:,j)),'s','color','#00539a','MarkerSize',mksz,'MarkerEdgeColor','#198038','MarkerFaceColor','#198038','linewidth',lw,'capsize',cpsz)
errorbar(x2,mean(we),mean(we)-min(we),max(we)-mean(we),'s','color','#00539a','MarkerSize',mksz,'MarkerEdgeColor','#a56eff','MarkerFaceColor','#a56eff','linewidth',lw,'capsize',cpsz)
hold off
set(gca, 'YGrid', 'on', 'XGrid', 'off')
ticklabel = get(gca,'TickLabel');
set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz2,'FontWeight','bold');
set(gca, 'XTick',[x1 x3 x2], 'XTickLabel',Cell{1})
xlim([x1-xgp x2+xgp])
title(dvd_title(j),'FontName','Consolas','fontsize',ftsz,'FontWeight','bold');

clc;
report_table = [median(dvd) median(we) median(pi_cf); mean(dvd) mean(we) mean(pi_cf);min(dvd) min(we) min(pi_cf);max(dvd) max(we) max(pi_cf)]';
fprintf('%25s %s\n','','**Debt valuation decomposition**')
fprintf('%25s %10s %10s %15s %15s\n','','median','mean','5th quantile','95th quantile')
for j = 1:6
fprintf('%25s %10.2f %10.2f %15.2f %15.2f\n',dvd_title(j),[report_table(j,:)])
end

%% Counterfactuals

% unwrap the structure
A_sampler = BVAR_est.A_sampler;
Sigma_sampler = BVAR_est.Sigma_sampler;
data_index = BVAR_est.data_index;
diff_index = BVAR_est.diff_index;

% time discount rate for Debt value decomposition
betat = .99;

%% counter-factuals
IRFcfQ = [];
%cfra = mean(we./dvd(:,4))-1;
cfra = .1;
CFRATIO = [1 1-cfra 1+cfra];
for jj = 1:length(CFRATIO)
% counterfactual sum of primary surplus path
cfratio = CFRATIO(jj);
ps_cf = cfratio*DVD(:,4);

% Now the PS follows an independent AR process
A_sampler_cf = A_sampler;
A_sampler_cf(:,ps,:) = 0;
A_sampler_cf(1+ps,ps,:) = betat;

% re-compute the impact
IRFcf = zeros(M,H+1,N);
impactcf = squeeze(IRF(:,1,:));
impactcf(Pi,:) = impactcf(Pi,:) + (DVD(:,4) - ps_cf)';
impactcf(ps,:) = ps_cf*(1-betat);

% compute the counterfactual IRF
for j = 1:N
    Acfj = A_sampler_cf(:,:,j);
    if regressor_constant
        Acfj = Acfj(2:end,:)';
    else
        Acfj = Acfj';
    end
    Acfj_comp = [Acfj; [kron(eye(p-1),eye(M)) zeros(M*(p-1),M)]];
    impactcf_comp = [impactcf(:,j); zeros(M*(p-1),1)];
    irfcfj = impactcf_comp;    
    for h = 0:H
        IRFcf(:,h+1,j) = irfcfj(1:M,:);
        irfcfj = Acfj_comp*irfcfj;
    end    
end
IRFcfQ = cat(3,IRFcfQ,median(IRFcf,3));
end

% variable name
y_nm_full = ["output";"consumption";"investment";"wage";"labor";"market value of debt";"primary surplus";"govt spending";"transfer";"inflation";"interest rate";"nominal bond return";"output gap";"primary surplus"];
y_nm = y_nm_full(data_index);
hrz = linspace(0,H,H+1)';
mycol = {'#00539a' '#ff832b' '#198038' '#a56eff'};
lw1 = 2.5;
lw2 = 1.5;
aph = 0.5;
fsz_tl = 24;
fsz_stl = 36;
fsz_lg = 12;
fsz_lb = 24;
fsz_tk = 24;
%%
lw1 = 4;
figure('name','Counterfactual fiscal backing','color','w','WindowState','maximized')
tiledlayout(1,3, 'Padding', 'none', 'TileSpacing', 'compact');    
for j = [7 5 1]
    nexttile
    hold on
    yline(0,'-','color',[1 1 1]*.25,'linewidth',1.5)
    plot(hrz,median(IRF(j,:,:),3)','color',mycol{1},'linewidth',lw1);
    plot(hrz,squeeze(IRFcfQ(j,:,1)),'--','color',mycol{2},'linewidth',lw1);
    plot(hrz,squeeze(IRFcfQ(j,:,2)),'-.','color',mycol{3},'linewidth',lw1);
    plot(hrz,squeeze(IRFcfQ(j,:,3)),':','color',mycol{4},'linewidth',lw1);
    hold off
    set(gca, 'YGrid', 'on', 'XGrid', 'off')
    ticklabel = get(gca,'TickLabel');
    set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',fsz_tk,'FontWeight','bold');
    title(y_nm(j),'FontName','Consolas','fontsize',fsz_tl,'FontWeight','bold');
    if j == 7
        ylabel('%')
        legend('',sprintf('%s','Estimated VAR'),sprintf('%s = %s','$\sum^\infty_{h=0}\frac{\partial s_{t+h}}{\partial \varepsilon_t}$','Estimated'),sprintf('%s = %.1f%s%s','$\sum^\infty_{h=0}\frac{\partial s_{t+h}}{\partial \varepsilon_t}$',1-cfra,'$\times$','Estimated'),sprintf('%s = %.1f%s%s','$\sum^\infty_{h=0}\frac{\partial s_{t+h}}{\partial \varepsilon_t}$',1+cfra,'$\times$','Estimated'),'interpreter','latex','location','south','FontSize',fsz_lg,'FontName','Consolas','FontWeight','bold');        
        legend boxoff
    end
end