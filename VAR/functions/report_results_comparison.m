function report_results_comparison(IRF,IRF_DSGE,data_index,conf_lev)

% measure sizes
[M H1 N] = size(IRF);
H_DSGE = size(IRF_DSGE,2)-1;

% horizon of IRF
H = H1-1;

% horizon to show
H_show = H;

% significance level
sig_lev = (1 - conf_lev) / 2;

% horizon for x-axis
hrz = linspace(0,H,H+1)';
hrz_dsge = linspace(0,H_DSGE,H_DSGE+1)';

% line width for plot
lw = 1.5;
lw2 = 2.5;
ftsz1 = 18;
ftsz2 = 16;

% variable name
y_nm_full = ["output";"consumption";"investment";"wage";"labor";"market value of debt";"primary surplus";"govt spending";"transfer";"inflation";"interest rate";"nominal bond return";"output gap";"primary surplus"];
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
    %p1 = plot(hrz,median(IRFj,2),'o-','color',mycol,'linewidth',lw,'markersize',5);    
    p2 = plot(hrz_dsge,IRF_DSGE(j,:)',':','color',mycol2,'linewidth',lw2,'markersize',5);
    hold off
    set(gca, 'YGrid', 'on', 'XGrid', 'off')
    ticklabel = get(gca,'TickLabel');
    set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz2,'FontWeight','bold');
    xticks([0:8:min(H,H_show)])
    xlim([0 min(H,H_show)])
    title(y_nm(j),'FontSize',ftsz1,'FontName','Consolas','FontWeight','bold');
    if j == 1
        legend([p1, p2],'VAR','DSGE','location','northwest','FontSize',ftsz1,'FontName','Consolas','FontWeight','bold');
        legend boxoff
end

end