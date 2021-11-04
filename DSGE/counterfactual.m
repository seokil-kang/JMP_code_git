% counterfactual exercise: changes transfer response to debt

% house keeping
clear;close all;clc;

% choose the version of model
cd benchmark

% load posterior sampler
load result/SMC_posterior
theta = squeeze(theta(:,:,end));
clearvars -except theta

% length of horizon for impulse response function
H = 40;

% confidence & significance level
conf_lev = .9;
sig_lev = (1-conf_lev)/2;

% add path for submodules
addpath('../functions','../gensys','model','../prior')

% prior info setup
[prior_info nm]= prior_info_setup;

% counterfactual setup

% which parameter to change?
cfparam_index = find(nm=="psi_b");

% what is the mean of the estimate?
param_bar = mean(theta(cfparam_index,:));

% how much to change?
cfparam1 = 5;
cfparam2 = 1/5;

% estimated parameter and counterfactuals
theta1 = theta;
theta1(cfparam_index,:) = cfparam1*theta1(cfparam_index,:);
theta2 = theta;        
theta2(cfparam_index,:) = cfparam2*theta2(cfparam_index,:);

% computation takes time a bit
if exist('result/counterfactuals.mat')
    load result/counterfactuals.mat
    
    % call variable indice
    variable_listing
    variable_listing_name
else
    tic;
    % compute IRF and Debt valuation decomposition
    [IRF0 DVD0] = IRF_and_DVD(theta,H);
    [IRF1 DVD1] = IRF_and_DVD(theta1,H);
    [IRF2 DVD2] = IRF_and_DVD(theta2,H);
    
    % call variable indice
    variable_listing
    variable_listing_name
    
    % select variables of interest
    var_select = [y infl Q mktb lambda R RL s];
    IRF0 = IRF0(var_select,:,:);
    IRF1 = IRF1(var_select,:,:);
    IRF2 = IRF2(var_select,:,:);
    
    % save the result
    save('result/counterfactuals.mat','IRF0', 'IRF1', 'IRF2', 'DVD0', 'DVD1', 'DVD2', 'cfparam_index', 'cfparam1', 'cfparam2','var_select')
    toc;
end

% clear up the submodule path
rmpath('../functions','../gensys','model','../prior')
cd ..
%% report result
mycol = {'#00539a' '#ff832b' '#198038' '#a56eff'};
linsty = {'-' ':' '--' '-.'};
mksty = {'none' 'none' 'none' 'none'};
%mksty = {'x' 'o' 'd' 's'};
lw1 = 4;
lw2 = 1.5;
ftsz = 20;
ftsz2 = 16;
apha = .5;

% x-axis horizon
hrz = linspace(0,H,H+1)';

% IRF figure
figure('name','DSGE Counterfactual IRF','color','w','WindowState','maximized')
tiledlayout(2,4, 'Padding', 'none', 'TileSpacing', 'compact');
for j = 1:8
    nexttile
    hold on
    yline(0,'color',[1 1 1]*apha,'linewidth',lw2)
    plot(hrz,squeeze(mean(IRF0(j,:,:),3)),'color',mycol{1},'linestyle',linsty{1},'marker',mksty{1},'linewidth',lw1)
    plot(hrz,squeeze(mean(IRF1(j,:,:),3)),'color',mycol{2},'linestyle',linsty{2},'marker',mksty{2},'linewidth',lw1)
    plot(hrz,squeeze(mean(IRF2(j,:,:),3)),'color',mycol{3},'linestyle',linsty{3},'marker',mksty{3},'linewidth',lw1)
    hold off
    set(gca, 'YGrid', 'on', 'XGrid', 'off')
    ticklabel = get(gca,'TickLabel');
    set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz,'FontWeight','bold');
    title(var_name_list(var_select(j)),'FontName','Consolas','fontsize',ftsz,'FontWeight','bold');
    if j == 1
        legend('',sprintf('%s%s = %.2f%s','\',sprintf('%s',nm(cfparam_index)),param_bar,'(estimated)'),sprintf('%s%s = %.2f','\',nm(cfparam_index),mean(theta1(cfparam_index,:))),sprintf('%s%s = %.2f','\',nm(cfparam_index),mean(theta2(cfparam_index,:))),'location','southeast','FontName','Consolas','fontsize',ftsz2,'FontWeight','bold');
        legend boxoff
    end
end

% figure debt valuation decomposition
figure('name','DSGE Counterfactual DVD','color','w','WindowState','maximized')
dvd_title  = ["Inflation at impact" "Bond price revaluation" "Discount rate effect" "Fiscal backing" "Wealth effect" "inflation coun-fact"];
tiledlayout(1,4, 'Padding', 'none', 'TileSpacing', 'compact');
for j = 1:4
    nexttile
    hold on
    [f0,x0] = ksdensity(rmoutliers(DVD0(:,j),'percentiles',[sig_lev 1-sig_lev]*1e2));
    [f1,x1] = ksdensity(rmoutliers(DVD1(:,j),'percentiles',[sig_lev 1-sig_lev]*1e2));
    [f2,x2] = ksdensity(rmoutliers(DVD2(:,j),'percentiles',[sig_lev 1-sig_lev]*1e2));
    p0 = plot(x0,f0,'color',mycol{1},'linestyle',linsty{1},'marker',mksty{1},'linewidth',lw1);
    p1 = plot(x1,f1,'color',mycol{2},'linestyle',linsty{2},'marker',mksty{2},'linewidth',lw1);
    p2 = plot(x2,f2,'color',mycol{3},'linestyle',linsty{3},'marker',mksty{3},'linewidth',lw1);
    ticklabel = get(gca,'TickLabel');
    set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz,'FontWeight','bold');
    title(dvd_title(j),'FontName','Consolas','fontsize',ftsz,'FontWeight','bold');
    if j == 1
        legend([p0 p1 p2],sprintf('%s%s = %.2f%s','\',sprintf('%s',nm(cfparam_index)),param_bar,'(estimated)'),sprintf('%s%s = %.2f','\',nm(cfparam_index),mean(theta1(cfparam_index,:))),sprintf('%s%s = %.2f','\',nm(cfparam_index),mean(theta2(cfparam_index,:))),'location','northwest','FontName','Consolas','fontsize',ftsz2,'FontWeight','bold');
        legend boxoff
    end
end