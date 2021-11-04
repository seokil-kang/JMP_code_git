% plot the IRF and Debt value decomposition to the contractionary MP shock

function report_counterfactual(IRF, DVD, conf_lev,cf)
% significance level
sig_lev = (1-conf_lev)/2;

% field names of structures
nm_irf = fieldnames(IRF);
nm_dvd = fieldnames(DVD);

% measure the sizes
m = length(nm_irf);
[N,nd] = size(DVD.(nm_dvd{1}));
[n_obs,H,~] = size(IRF.(nm_irf{1}));

% IRF horizon length
H = H-1;

% x-axis for IRF
hrz = linspace(0,H,H+1)';

% some labels
y_nm = ["Output" "Investment" "Wage" "Labor" "Gov't Debt" "Primary Surplus" "Inflation" "Interest Rate" "Holding Return"]';
d_nm = ["Inflation impact" "Bond price revaluation" "Sum of expected discount rate" "Sum of expected primary surplus"]';

% IRF
lw = 2.5;
figure('name','IRF to Contractionary MP Shock','color','w')
for j = 1:n_obs
    subplot(3,3,j)
    hold on
    yline(0,'r')
    for jj = 1:m
        irfqj = squeeze(quantile(IRF.(nm_irf{jj})(j,:,:),[.5 sig_lev 1-sig_lev],3));
        if jj == 1
            p0 = plot(hrz,irfqj(:,1),'d-','linewidth',lw);
        elseif jj == 2
            p1 = plot(hrz,irfqj(:,1),'x--','linewidth',lw);
        else
            p2 = plot(hrz,irfqj(:,1),'o-.','linewidth',lw);
        end
        %plot(hrz,irfqj(:,2:3),'--b','linewidth',1.5);        
    end
    hold off
    xlim([0 H])
    title(y_nm(j))
    set(gca,'YGrid','on')
end
legend([p0,p1,p2],sprintf('$$%s$$ = %.1f(%s)',cf.nm_cf,cf.cf0,'estimated'),sprintf('$$%s$$ = %.1f',cf.nm_cf,cf.cf1),sprintf('$$%s$$ = %.1f',cf.nm_cf,cf.cf2),'location','southeast','interpreter','latex')
legend boxoff
%%
% Debt valuation decomposition
figure('name','Debt Valuation Decomposition to Contractionary MP Shock','color','w')
for j = 1:nd
subplot(1,nd,j)
hold on
for jj = 1:m
    dj = rmoutliers(DVD.(nm_dvd{jj})(:,j),'percentiles',[sig_lev 1-sig_lev]*1e2);
    [fj xj] = ksdensity(dj);
    if jj == 1
        plot(xj,fj,'linewidth',lw)
    elseif jj == 2
        plot(xj,fj,'--','linewidth',lw)
    elseif jj == 3
        plot(xj,fj,'-.','linewidth',lw)
    end
    %histogram(dj,'normalization','pdf')
end
hold off
title(d_nm(j))
end
legend(sprintf('$$%s$$ = %.1f(%s)',cf.nm_cf,cf.cf0,'estimated'),sprintf('$$%s$$ = %.1f',cf.nm_cf,cf.cf1),sprintf('$$%s$$ = %.1f',cf.nm_cf,cf.cf2),'location','northeast','interpreter','latex')
legend boxoff
end