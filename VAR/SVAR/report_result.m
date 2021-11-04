clear; clc; close all;
load results/results
conf_lev = .9;
H = 40;
[k,M,N] = size(Aplustilde);
p = (k-1)/M;
IRF = zeros(M,H+1,N);
DVD = zeros(N,4);
FETV = zeros(M,H+1,N);
IDSV = zeros(M,H+1,N);
MP_rule = zeros(N,2);
cd ../functions
variable_indexing(data_index)
cd ../SVAR
if exist('y')
    output = y;
elseif exist('x')
    output = x;
end
for j = 1:N
    Ainvj = inv(A0tilde(:,:,j));
    Bj = Aplustilde(:,:,j)*Ainvj;
    Bcompj = [Bj(1:end-1,:)'; [eye(M*p-M) zeros(M*p-M,M)]];
    impactj = Ainvj(1,:)';
    %impactj = impactj / impactj(R) * .25;
    %impactj = impactj * 2.5;
    
    if impactj(R) < 0
        impactj = -impactj;
    end
    
    %{
    irfj = [];    
    fetv = Sigma_sampler(:,:,j);
    idsv = impactj*impactj';                
    for h = 0:H
        Bh = Bcompj^h;
        irfjh = Bh(1:M,1:M)*impactj;
        irfj = [irfj irfjh];
        FETV(:,h+1,j) = diag(Bh(1:M,1:M)*fetv*Bh(1:M,1:M)');
        IDSV(:,h+1,j) = diag(Bh(1:M,1:M)*idsv*Bh(1:M,1:M)');
    end
        %}        
    irfj = [];
    irfjh = [impactj;  zeros(M*(p-1),1)];
    fetv = kron(reshape(double(1:p*p == 1),p,p),Sigma_sampler(:,:,j));
    idsv = kron(reshape(double(1:p*p == 1),p,p),impactj*impactj');
    for h = 0:H
        irfj = [irfj irfjh(1:M)];
        irfjh = Bcompj * irfjh;
        FETV(:,h+1,j) = diag(fetv(1:M,1:M));    
        fetv = Bcompj * fetv * Bcompj';
        IDSV(:,h+1,j) = diag(idsv(1:M,1:M));    
        idsv = Bcompj * idsv * Bcompj';
    end
    IRF(:,:,j) = irfj;
    
    betat = .99;
    impactjvec = [impactj;  zeros(M*(p-1),1)];
    
    cPS = zeros(M*p,1);
    cPS(6) = 1;
    cPi = zeros(M*p,1);
    cPi(7) = 1;
    cRB = zeros(M*p,1);
    cRB(9) = 1;
    
    revaluation_effect = cRB'*impactjvec;
    inflation_distribution = cPi'*impactjvec;
    debt_value = revaluation_effect - inflation_distribution;
    expected_PS_path = cPS' * inv(eye(M*p) - betat*Bcompj)*impactjvec;
    expected_real_DC = -(cRB' - cPi') * inv(eye(M*p) - betat*Bcompj)*impactjvec + debt_value;
    
    DVD(j,:) = [inflation_distribution revaluation_effect expected_real_DC expected_PS_path];
    
    % MP rule coefficients
    
    MP_rule(j,:) = -[A0tilde(output,1,j) A0tilde(Pi,1,j)]/A0tilde(R,1,j);
end
% [Ainv(1,8,j) Ltilde(1,8,1,j)] are the same for all j
% they are the impact response of R on MP shock

% cumulate the IRF for differenced data
if length(diff_index) > 0
    IRF(diff_index,:,:) = cumsum(IRF(diff_index,:,:),2);
end

% Forecast Error decomposition
FEVD = cumsum(IDSV,2) ./ cumsum(FETV,2);

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

%%
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