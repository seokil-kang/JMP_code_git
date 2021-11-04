function report_results(IRF,DVD,FEVD,data_index,conf_lev)

% measure sizes
[M H1 N] = size(IRF);

% horizon of IRF
H = H1-1;

% horizon to show
H_show = 40;

% significance level
sig_lev = (1 - conf_lev) / 2;

% horizon for x-axis
hrz = linspace(0,H,H+1)';

% line width for plot
lw = 1.5;
ftsz1 = 18;
ftsz2 = 16;

% variable name
y_nm_full = ["output";"consumption";"investment";"wage";"labor";"debt";"surplus";"govt spending";"transfer";"inflation";"interest rate";"holding return";"output gap";"surplus"];
y_nm = y_nm_full(data_index);

%{
figure('name','IRF: Sign Restriction','color','w')
for j = 1:M
    IRFj = squeeze(IRF(j,:,:));
    IRFjQ = quantile(IRFj,[.5 sig_lev 1-sig_lev],2);
    subplot(3,3,j)
    hold on
    yline(0,'r')
    plot(hrz,IRFjQ(:,1),'.-k','linewidth',lw)
    plot(hrz,IRFjQ(:,2:3),':b','linewidth',lw)
    hold off
    xlim([0 H])
    title(y_nm(j))
end
%}

N_quantile = 2;
posterior_band_q = linspace(sig_lev/2,.5,N_quantile);
%posterior_band_q = [sig_lev/2 1/6 .5];
%N_quantile = 3;
mycol = [0 .45 .74];
mycol = [.0 .0 .5];
figure('WindowState','maximized','name','IRF to MP Shock','color','w')
tiledlayout(3,3, 'Padding', 'none', 'TileSpacing', 'compact');
for j = 1:M
    IRFj = squeeze(IRF(j,:,:));
    nexttile
    hold on
    yline(0,'r')
    for jj = 1:N_quantile
        low_quantile = quantile(IRFj,1-posterior_band_q(jj),2);
        high_quantile = quantile(IRFj,posterior_band_q(jj),2);
        if jj == 1
            fill([hrz', flip(hrz')], [low_quantile' flip(high_quantile')],mycol,'facealpha',1/4,'linestyle','none');
        else
            fill([hrz', flip(hrz')], [low_quantile' flip(high_quantile')],mycol,'facealpha',jj/(N_quantile^1.2),'linestyle','none');
        end
    end
    %plot(hrz,mean(IRFj,2),'.-','color',mycol,'linewidth',.1,'markersize',10);
    plot(hrz,mean(IRFj,2),'o-','color',mycol,'linewidth',lw,'markersize',5);
    hold off
    set(gca, 'YGrid', 'on', 'XGrid', 'off')
    ticklabel = get(gca,'TickLabel');
    set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz2,'FontWeight','bold');
    xticks([0:8:min(H,H_show)])
    xlim([0 min(H,H_show)])
    title(y_nm(j),'FontSize',ftsz1,'FontName','Consolas','FontWeight','bold');
end

d_nm = ["Inflation impact";"Bond price revaluation";"Sum of expected discount rate";"Sum of expected primary surplus"];
figure('WindowState','maximized','name','Debt Value Decomposition to MP Shock','color','w')
tiledlayout(1,4, 'Padding', 'none', 'TileSpacing', 'compact');
for j = 1:4
    Dj = rmoutliers(DVD(:,j),'percentiles',[sig_lev 1-sig_lev]*100);
    nexttile
    hold on
    histogram(Dj,'normalization','pdf','facecolor',mycol,'facealpha',3/4);%,'linestyle','none')
    %[fj xj] = ksdensity(Dj);
    %plot(xj,fj,'linewidth',3,'color',mycol)
    %fill([xj' flip(xj')], [fj' flip(fj')],mycol,'facealpha',2/4,'linestyle','none')
    hold off
    ticklabel = get(gca,'TickLabel');
    set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz2,'FontWeight','bold');
    title(d_nm(j),'FontSize',ftsz1,'FontName','Consolas','FontWeight','bold');
end
xlabel('% change','FontSize',18,'FontName','Consolas');

colr = [mycol
    0.50 0.50 0.50
 0.64 0.08 0.18
 1.00 0.41 0.16
 0.93 0.69 0.13
 0.47 0.67 0.19
 0.30 0.75 0.93
 0.08 0.18 0.49
 0.49 0.18 0.56];

% Forecast Error Decomposition
figure('WindowState','maximized','name','Forecast Error Decomposition: MP Shock','color','w')
tiledlayout(3,3, 'Padding', 'none', 'TileSpacing', 'compact');
for j = 1:M
    nexttile
    bar(hrz,squeeze(mean(FEVD(j,:,:),3)),1,'stacked','facealpha',3/4)
    colororder(colr)
    title(y_nm(j))
    set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz2,'FontWeight','bold');
    xticks([0:8:min(H,H_show)])
    xlim([0 min(H,H_show)])
    title(y_nm(j),'FontSize',ftsz1,'FontName','Consolas','FontWeight','bold');
end
end