% report parameter estimation of SMC

% house keeping
clear;close all;clc;

% confidence level
conf_lev = .9;

% choose the version of model
cd benchmark

% load posterior draws
load result/SMC_posterior

% measure the posterior sampler
[n N N_stage] = size(theta);

% add path for submodules
addpath('../functions','../gensys','model','../prior')

% prior info setup
[prior_info nm nmtx]= prior_info_setup;
cd ..

% prior sampler
theta_prior = theta(:,:,1);

% posterior sampler
theta_posterior = theta(:,:,end);

% two-sided significance level
sig_lev = (1-conf_lev)/2;

%% report results
mycol = {'#00539a' '#ff832b' '#198038' '#a56eff'};
lw = 3;
ftsz = 24;
ftsz2 = 20;

%%
figure('name','SMC hyperparameter adaptation','color','w','WindowState','maximized')
tiledlayout(3,1, 'Padding', 'none', 'TileSpacing', 'compact');
nexttile
hold on
yline(.25,'--')
plot(Acceptance_rate,'color',mycol{1},'linewidth',lw)
hold off
set(gca, 'YGrid', 'on', 'XGrid', 'off')
ticklabel = get(gca,'TickLabel');
set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz2,'FontWeight','bold');
xlim([1 500])
title('Mutation Acceptance Rate','FontName','Consolas','fontsize',ftsz,'FontWeight','bold');
nexttile
plot(c,'color',mycol{1},'linewidth',lw)
set(gca, 'YGrid', 'on', 'XGrid', 'off')
ticklabel = get(gca,'TickLabel');
set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz2,'FontWeight','bold');
xlim([0 500])
title('Scale Parameter for Mutation','FontName','Consolas','fontsize',ftsz,'FontWeight','bold');
nexttile
plot(ESS,'color',mycol{1},'linewidth',lw)
set(gca, 'YGrid', 'on', 'XGrid', 'off')
ticklabel = get(gca,'TickLabel');
set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz2,'FontWeight','bold');
xlim([1 500])
title('Effective Sample Size','FontName','Consolas','fontsize',ftsz,'FontWeight','bold');
xlabel('stage','FontName','Consolas','fontsize',ftsz,'FontWeight','bold');

%% plot sampler density
figure('name','SMC Sampler density','color','w','WindowState','maximized')
if n < 37
    tiledlayout(6,6, 'Padding', 'none', 'TileSpacing', 'compact');
else
    tiledlayout(6,7, 'Padding', 'none', 'TileSpacing', 'compact');
end
for j = 1:n
    nexttile
    [f1 x1] = ksdensity(theta_prior(j,:));
    [f2 x2] = ksdensity(theta_posterior(j,:));    
    hold on
    plot(x1,f1,'color',mycol{1},'linewidth',lw);
    plot(x2,f2,'color',mycol{2},'linewidth',lw);
    hold off
    ticklabel = get(gca,'TickLabel');
    set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz2,'FontWeight','bold');
    title(nmtx(j),'FontName','Consolas','fontsize',ftsz,'FontWeight','bold','interpreter','latex');
end

%%
ysct = theta_posterior(find(nm=='phi_pi'),:);
corr_param = corr(ysct',theta_posterior')';
[~, corr_index] = sort(abs(corr_param),'descend');
figure('name','Posterior Sampler Correlation','color','w','WindowState','maximized')
if n < 37
    tiledlayout(6,6, 'Padding', 'none', 'TileSpacing', 'compact');
else
    tiledlayout(6,7, 'Padding', 'none', 'TileSpacing', 'compact');
end
for j = corr_index'
    nexttile
    xsct = theta_posterior(j,:);        
    p0 = scatter(xsct,ysct,'.','MarkerEdgeColor',mycol{1});
    coef = polyfit(xsct,ysct,1);
    hline = refline(coef(1),coef(2));    
    hline.Color = [255 131 43]/255*abs(corr(xsct',ysct'));
    set(hline,'linewidth',2)
    ticklabel = get(gca,'TickLabel');
    set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz2,'FontWeight','bold');
    title(sprintf('%s(%.2f)',nmtx(j),corr_param(j)),'FontName','Consolas','fontsize',ftsz,'FontWeight','bold','interpreter','latex');
end

%%
stages = round(linspace(1,N_stage,50));
n_stage = length(stages);
figure('name','Particle Sampler Histogram','color','w','WindowState','maximized')
if n < 37
    tiledlayout(6,6, 'Padding', 'none', 'TileSpacing', 'compact');
else
    tiledlayout(6,7, 'Padding', 'none', 'TileSpacing', 'compact');
end
for j = 1:n
    nexttile
    xmat = [];
    zmat = [];
    for each_stage = stages
        [fn xn] = ksdensity(theta(j,:,each_stage));
        xmat = [xmat;xn];
        zmat = [zmat;fn];
    end
    xmat = xmat';
    zmat = zmat';
    ymat = repmat(stages,numel(xn),1);
    mesh(xmat,ymat,zmat,'facecolor','flat');
    colormap turbo
    ticklabel = get(gca,'TickLabel');
    set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz2,'FontWeight','bold');
    title(nmtx(j),'FontName','Consolas','fontsize',ftsz,'FontWeight','bold','interpreter','latex');
    ylim([0 size(theta,3)])
end

%%
clc;
% report the credible interval
report_table = [nm prior_info(:,[2 3]) quantile(theta_prior,[sig_lev 1-sig_lev],2) mean(theta_posterior,2) std(theta_posterior')' quantile(theta_posterior,[sig_lev 1-sig_lev],2)];
fprintf('%56s\n','SMC estimation result table')
fprintf('%32s %32s\n','prior', 'posterior')
fprintf('%15s %5s %6s %10s %14s %6s %10s\n','parameter','mean', 'stdev', '90 CS','mean', 'stdev', '90 CS')
fprintf('%15s %5.2f %6.2f [%5.2f, %5.2f] %10.2f %6.2f [%5.2f, %5.2f]\n',report_table')

% compute the marginal data density
log_marginal_data_density = sum(log(diag(w_tilde' * W(:,1:end-1))/N));
fprintf('log marginal data denstiy = %.0f\n',log_marginal_data_density)