% house keeping
clear;clc;close all;

% addpath for required functions
addpath('../functions')

% Length of IRF horizon
H = 40;

% confidence level
conf_lev = .9;

% load the BPSVAR estimates
load results/BVAR_posterior

% call variable indexing
variable_indexing(BVAR_est.data_index);

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

% denote the policy index
policy_index = R;

% IRF and Debt value decomposition basket
IRF = zeros(M,H+1,N);
DVD = zeros(N,4);
FETV = zeros(M,H+1,N);
IDSV = zeros(M,H+1,N);

% time discount rate for Debt value decomposition
betat = .99;

% basis vectors for variable selection
cps = double(1:M*p == ps)';
cpi = double(1:M*p == Pi)';
cRB = double(1:M*p == RB)';
cb = double(1:M*p == b)';

% draw IRF
for j = 1:N
    
    % select the coefficient estimates
    Aj = BVAR_est.A_sampler(:,:,j);
    
    % select the identified shock impact
    impactj = BVAR_est.impact(:,j);
        
    if impactj(policy_index) < 0
        impactj = -impactj;
    end
    
    % normalize to 25 basis-points shock
    impactj = impactj / impactj(policy_index) * .25;
    
    % is there a constant vector in VAR?
    if regressor_constant
        Aj = Aj(1:end-1,:)';
    else
        Aj = Aj';
    end
    
    % companion VAR
    Aj_comp = [Aj; [kron(eye(p-1),eye(M)) zeros(M*(p-1),M)]];
    
    % compute the IRF beyond the restricted horizon
    irfj = [];
    fetv = BVAR_est.Sigma_sampler(:,:,j);
    idsv = impactj*impactj';
    for h = 0:H
        Ah = Aj_comp^h;
        irfj = [irfj Ah(1:M,1:M)*impactj];
        FETV(:,h+1,j) = diag(Ah(1:M,1:M)*fetv*Ah(1:M,1:M)');
        IDSV(:,h+1,j) = diag(Ah(1:M,1:M)*idsv*Ah(1:M,1:M)');
    end
    IRF(:,:,j) = irfj;
    
    % impact vector for debt value decomposition
    impactj_vec = [impactj; zeros(M*(p-1),1)];
    
    % debt value decomposition
    inflation_impact = cpi' * impactj_vec;
    bond_price_reval = cRB' * impactj_vec;
    debt_reval = bond_price_reval - inflation_impact;
    Ex_dc_path = -(cRB - cpi)' * inv(eye(M*p) - betat * Aj_comp) * impactj_vec + debt_reval;
    Ex_ps_path = cps' * inv(eye(M*p) - betat * Aj_comp) * impactj_vec;
    
    % collect the decomposition
    DVD(j,:) = [inflation_impact bond_price_reval Ex_dc_path Ex_ps_path];
end

% cumulate the IRF for differenced data
if length(BVAR_est.diff_index) > 0
    IRF(BVAR_est.diff_index,:,:) = cumsum(IRF(BVAR_est.diff_index,:,:),2);
end

% Forecast Error decomposition
FEVD = cumsum(IDSV,2) ./ cumsum(FETV,2);

%%
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
y_nm = y_nm_full(BVAR_est.data_index);

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


dvd_title  = ["Inflation at impact" "Bond price revaluation" "Discount rate effect" "Fiscal backing" "Wealth effect" "inflation coun-fact"];
dvd = rmoutliers(DVD,'percentiles',[(1-conf_lev)/2 1-(1-conf_lev)/2]*1e2,1);
we = dvd(:,2) - dvd(:,1) - dvd(:,3);
pi_cf = dvd(:,1) - (we-dvd(:,4));
report_table = [median(dvd) median(we) median(pi_cf); mean(dvd) mean(we) mean(pi_cf);min(dvd) min(we) min(pi_cf);max(dvd) max(we) max(pi_cf)]';
fprintf('%25s %s\n','','**Debt valuation decomposition**')
fprintf('%25s %10s %10s %15s %15s\n','','median','mean','5th quantile','95th quantile')
for j = 1:6
fprintf('%25s %10.2f %10.2f %15.2f %15.2f\n',dvd_title(j),[report_table(j,:)])
end