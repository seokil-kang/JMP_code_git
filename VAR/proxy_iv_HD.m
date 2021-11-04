% Estimate IRF with proxy instrument variable

% house keeping
clear;clc;close all;

% addpath for required functions
addpath('functions')

% confidence level
conf_lev = .9;

% significance level
sig_lev = (1-conf_lev)/2;

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

% measure the sample size for IV estimation
T = length(Z);

% unwrap the structure
A_sampler = BVAR_est.A_sampler;
Sigma_sampler = BVAR_est.Sigma_sampler;
data_index = BVAR_est.data_index;
diff_index = BVAR_est.diff_index;

if exist('HD_VAR_benchmark.mat')
    load HD_VAR_benchmark    
else    
    % basket for relevance of proxy IV & MP Shock
    Phi = [];
    MP_shock_hat = [];
    Y_HD = [];
    
    for j = 1:N
        
        % retrieve coefficient matrix from vectorization
        Aj = A_sampler(:,:,j);
        
        % residual of reduced form VAR
        Uj = Y - X * Aj;
        
        % draw Sigma from sampler
        Sigmaj = Sigma_sampler(:,:,j);
        
        % recover the IV correlation
        phij = sqrt((Z'*Uj/T)*inv(Uj'*Uj)*(Uj'*Z));    
        
        % estimate the impact
        impactj = (Uj'*Z/T)/phij;
        
        % recover the monetary policy shock
        MP_shock_hatj = (impactj'*inv(Sigmaj)*Uj' / phij)';
        
        % is there a constant vector in VAR?
        if regressor_constant
            Aj = Aj(2:end,:)';
        else
            Aj = Aj';
        end
        
        % companion VAR
        Aj_comp = [Aj; [kron(eye(p-1),eye(M)) zeros(M*(p-1),M)]];
        
        % companion impact vector
        impact_compj = [impactj; zeros(M*(p-1),1)];
        
        % Historical decomposition   
        HDj = [];
        y_mp = [];
        for h = 0:T-2
            hdj = Aj_comp^h*impact_compj;
            HDj = [HDj hdj(1:M)];
            y_mpt = sum(HDj .* flip(MP_shock_hatj(2:h+2))',2);
            y_mp = [y_mp y_mpt];
        end
        
        % collect the estimates
        Phi = [Phi; phij];
        MP_shock_hat = [MP_shock_hat MP_shock_hatj];
        Y_HD = cat(3,Y_HD,y_mp);
    end
end
%%
close all;
y_nm_full = ["output";"consumption";"investment";"wage";"labor";"market value of debt";"primary surplus";"govt spending";"transfer";"inflation";"interest rate";"nominal bond return";"output gap";"primary surplus"];
y_nm = y_nm_full(data_index);
PS = Y(2:end,6);
PS_HD = squeeze(Y_HD(6,:,:));
PS_HDQ = quantile(PS_HD,[sig_lev 1-sig_lev .5],2);
new_measure = PS - PS_HD;
new_measureQs = quantile(new_measure,[sig_lev 1-sig_lev .5],2);
new_measureQ = new_measureQs(:,3);
lw = 2.5;
lw2 = 2;
ftsz = 24;
mycol = {'#00539a' '#ff832b' '#198038' '#a56eff' '#f1c21b' '#ff8389'};
mycols = [0 83 154; 255 131 43; 25 128 56; 165 110 255; 170 170 0]/255;
linsty = {'-' ':' '-.' '-.'};
mksty = {'none' 'none' 'none' 'none'};
%mksty = {'x' 'o' 'd' 's'};
figure('WindowState','maximized','name','Primary Surplus Without Fiscal Backing','color','w')
tiledlayout(1,1, 'Padding', 'none', 'TileSpacing', 'compact');
nexttile
hold on
%fill([dateq(2:end)', flip(dateq(2:end)')], [PS_HDQ(:,1)' flip(PS_HDQ(:,2)')],mycols(1,:),'facealpha',1/10,'linestyle','none');
p3 = bar(dateq(2:end),PS_HDQ(:,3),1,'edgecolor','none','facecolor',mycols(5,:),'facealpha',3/4);
p1 = plot(dateq(2:end),PS,'color',mycol{1},'linestyle',linsty{1},'marker',mksty{1},'linewidth',lw);
p2 = plot(dateq(2:end),new_measureQ,'color',mycol{6},'linestyle',linsty{2},'marker',mksty{2},'linewidth',3);
yline(0,':','color',mycols(3,:),'linewidth',1)    
hold off
recessionplot
axis tight
set(gca, 'YGrid', 'on', 'XGrid', 'off')    
ticklabel = get(gca,'TickLabel');
set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz,'FontWeight','bold');
legend([p1, p2 p3],'data','filtered','Fiscal backing','location','northwest','FontSize',16,'FontName','Consolas','FontWeight','bold');
legend boxoff
%%
nlag = 10;
mksz = 7;
figure('WindowState','maximized','name','Cross Correlation: Primary Surplus Without Fiscal Backing','color','w')
tiledlayout(3,3, 'Padding', 'none', 'TileSpacing', 'compact');
for j = 1:M
    nexttile
    [xcf1,lags,bounds1]=crosscorr(PS,Y(2:end,j),'numlags',nlag);
    [xcf2,lags,bounds2]=crosscorr(new_measureQ,Y(2:end,j),'numlags',nlag);
    hold on
    %yline(bounds1(1),':','color',mycols(3,:),'linewidth',2)    
    stem(lags,xcf1,'filled','-x','color',mycols(1,:),'markersize',mksz+5,'linewidth',lw2)
    stem(lags,xcf2,'filled','--s','color',mycols(2,:),'markersize',mksz,'linewidth',lw2)
    hold off
    set(gca, 'YGrid', 'on', 'XGrid', 'off')    
    ticklabel = get(gca,'TickLabel');
    set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz,'FontWeight','bold');    
    title(y_nm(j),'FontName','Consolas','fontsize',ftsz,'FontWeight','bold');
end
legend('original','filtered','location','northeast','FontSize',16,'FontName','Consolas','FontWeight','bold');
legend boxoff
%%
prxmty = Phi./std(MP_shock_hat)'/std(Z);
figure('WindowState','maximized','name','Proximity of Instrumental Variable','color','w')
tiledlayout(1,1, 'Padding', 'none', 'TileSpacing', 'compact');
nexttile
hold on
%histogram(rmoutliers(prxmty,'percentiles',[.5 99.5]),'normalization','pdf')
histogram(Phi,'normalization','pdf','edgecolor',mycols(1,:),'facecolor',mycols(1,:),'facealpha',.8)
histogram(prxmty,'normalization','pdf','edgecolor',mycols(2,:),'facecolor',mycols(2,:),'facealpha',.8)
hold off
ticklabel = get(gca,'TickLabel');
set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz,'FontWeight','bold');
legend('$\hat{\phi}$','$\frac{\hat{\phi}}{s(Z_t) s(\hat{\varepsilon^M_t})}$','interpreter','latex','FontName','Consolas','fontsize',60,'FontWeight','bold','location','north');
legend boxoff
xlim([.35 1])
if exist('HD_VAR_benchmark.mat')
elseif exist('HD_VAR_accounting.mat')
else
    save('HD_VAR_benchmark.mat','Phi','MP_shock_hat','Y_HD')    
end