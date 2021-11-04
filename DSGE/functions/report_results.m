function report_results(IRF,DVD,FEVD,var_select,conf_lev)

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
variable_listing_name;
e_nm = ["MP" "growth" "preference" "investment" "wage" "price" "spending" "transfer"];

N_quantile = 2;
posterior_band_q = linspace(sig_lev/2,.5,N_quantile);
mycol = [0 .45 .74];
mycol = [.0 .0 .5];
mycol = [0, 83, 154]/255/1.5;
figure('WindowState','maximized','name','IRF to MP Shock','color','w')
if M == 9        
        tiledlayout(3,3, 'Padding', 'none', 'TileSpacing', 'compact');
    else
        %tiledlayout(1,M, 'Padding', 'none', 'TileSpacing', 'compact'); 
        tiledlayout(2,ceil(M/2), 'Padding', 'none', 'TileSpacing', 'compact');
        %tiledlayout(ceil(M/2),2, 'Padding', 'none', 'TileSpacing', 'compact'); 
    end    
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
    plot(hrz,mean(IRFj,2),'o-','color',mycol,'linewidth',lw,'markersize',5);
    hold off    
    set(gca, 'YGrid', 'on', 'XGrid', 'off')
    ticklabel = get(gca,'TickLabel');
    set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz2,'FontWeight','bold');
    xticks([0:8:min(H,H_show)])
    xlim([0 min(H,H_show)])
    title(var_name_list{var_select(j)},'FontSize',ftsz1,'FontName','Consolas','FontWeight','bold');
end

d_nm = ["Inflation impact";"Bond price revaluation";"Sum of expected discount rate";"Sum of expected primary surplus"];
figure('WindowState','maximized','name','Debt Value Decomposition to MP Shock','color','w')
tiledlayout(1,4, 'Padding', 'none', 'TileSpacing', 'compact');
dvd = rmoutliers(DVD,'percentiles',[sig_lev 1-sig_lev]*100,1);
meanDj = mean(dvd);
minDj = min(dvd);
maxDj = max(dvd);
medianDj = median(dvd);
for j = 1:4
    nexttile
    hold on
    histogram(dvd(:,j),'normalization','pdf','facecolor',mycol)
    %[fj xj] = ksdensity(dvd(:,j));
    %plot(xj,fj,'linewidth',3,'color',mycol)
    %fill([xj' flip(xj')], [fj' flip(fj')],mycol,'facealpha',1/4,'linestyle','none')
    hold off
    set(gca,'LooseInset',get(gca,'TightInset'))
    title(d_nm(j),'FontSize',ftsz1,'FontName','Consolas','FontWeight','bold');
end
xlabel('% change','FontSize',ftsz1,'FontName','Consolas','FontWeight','bold');

report_table = [medianDj' meanDj' minDj' maxDj'];
fprintf('%35s %s\n','','**Debt valuation decomposition**')
fprintf('%35s %10s %10s %15s %15s\n','','median','mean','5th quantile','95th quantile')
for j = 1:4
fprintf('%35s %10.2f %10.2f %15.2f %15.2f\n',d_nm(j),[report_table(j,:)])
end

colr = [0.50 0.50 0.50
 0.64 0.08 0.18
 1.00 0.41 0.16
 0.93 0.69 0.13
 0.47 0.67 0.19
 0.30 0.75 0.93
 0.08 0.18 0.49
 0.49 0.18 0.56];

% Forecast Error Decomposition
figure('WindowState','maximized','name','Forecast Error Decomposition: MP Shock','color','w')
if M == 9        
        tiledlayout(3,3, 'Padding', 'none', 'TileSpacing', 'compact');
    else
        %tiledlayout(1,M, 'Padding', 'none', 'TileSpacing', 'compact'); 
        tiledlayout(2,ceil(M/2), 'Padding', 'none', 'TileSpacing', 'compact');
        %tiledlayout(ceil(M/2),2, 'Padding', 'none', 'TileSpacing', 'compact'); 
    end    
for j = 1:M
    nexttile
    bar(hrz(1:min(H,40)+1),squeeze(mean(FEVD(j,:,[end 1:end-1],:),4)),1,'stacked')
    colororder(colr)
    set(gca,'LooseInset',get(gca,'TightInset'))
    ticklabel = get(gca,'TickLabel');
    set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz1,'FontWeight','bold');
    xlim([0 min(H,40)])
    ylim([0 1])
    title(var_name_list{var_select(j)},'FontSize',ftsz1,'FontName','Consolas','FontWeight','bold');
end
legend(e_nm)
end