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
load ../DSGE/benchmark/result/Variance_Decomposition
FEVD_DSGE = squeeze(FEVD(:,:,end,:));

% unwrap the structure
A_sampler = BVAR_est.A_sampler;
Sigma_sampler = BVAR_est.Sigma_sampler;
data_index = BVAR_est.data_index;
diff_index = BVAR_est.diff_index;

% call the variable indexing
variable_indexing(data_index);

% measure the degree of freedom and # of draws
[K,~,N] = size(A_sampler);

% is there a constant?
if mod(K,M) == 0
    regressor_constant = 0;
elseif mod(K,M) == 1
    regressor_constant = 1;
end 

% lag length is...
p = (K-regressor_constant) / M;

% Variance decomposition basket
FETV = zeros(M,H+2,N);
IDSV = zeros(M,H+2,N);

% non policy indice
q_index = setdiff(1:M,policy_index);

% draw IRF
for j = 1:N
    % retrieve coefficient matrix from vectorization
    Aj = A_sampler(:,:,j);
    
    % residual of reduced form VAR
    Uj = Y - X * Aj;
    
    % policy residuals and non-policy residuals
    upj = Uj(:,policy_index);
    uqj = Uj(:,q_index);
    
    % draw Sigma from sampler
    Sigmaj = Sigma_sampler(:,:,j);
    
    % reordered Sigma
    Sigmaj_ro = Sigmaj([policy_index q_index],[policy_index q_index]);       
    
    % measure the sample size for IV estimation
    T = length(Z);
    
    % recover the IV correlation
    phij = sqrt((Z'*Uj/T)*inv(Uj'*Uj)*(Uj'*Z));
    
    % estimate the impact
    impact = (Uj'*Z/T)/phij;
    
    % contractionary shock
    if impact(policy_index) < 0
        impact = -impact;
    end
    
    % is there a constant vector in VAR?
    if regressor_constant
        Aj = Aj(2:end,:)';
    else
        Aj = Aj';
    end
    
    % companion VAR
    Aj_comp = [Aj; [kron(eye(p-1),eye(M)) zeros(M*(p-1),M)]];
    
    % companion impact vectors
    impact_comp = [impact; zeros(M*(p-1),1)];
    
    % compute the IRF beyond the restricted horizon    
    forcast_error_total_variance = kron(reshape(double(1:p*p == 1),p,p),Sigmaj);
    ID_shock_variance = impact_comp * impact_comp';
    for h = 0:H        
        FETV(:,h+1,j) = diag(forcast_error_total_variance(1:M,1:M));
        IDSV(:,h+1,j) = diag(ID_shock_variance(1:M,1:M));
        forcast_error_total_variance = Aj_comp*forcast_error_total_variance*Aj_comp';
        ID_shock_variance = Aj_comp * ID_shock_variance * Aj_comp';
    end
    LRFEVD = inv(eye(K-1) - Aj_comp) * kron(reshape(double(1:p*p == 1),p,p),Sigmaj) * inv(eye(K-1) - Aj_comp)';
    LRIDSV = inv(eye(K-1) - Aj_comp) * impact_comp * impact_comp' * inv(eye(K-1) - Aj_comp)';
    FETV(:,end,j) = diag(LRFEVD(1:M,1:M));
    IDSV(:,end,j) = diag(LRIDSV(1:M,1:M));
    
end

% Forecast Error decomposition
FEVD = cumsum(IDSV,2) ./ cumsum(FETV,2);
FEVD(:,end,:) = IDSV(:,end,:) ./ FETV(:,end,:);
%%
vd_interval = [[1:8:H]';H+1;H+2];
var_select = [y Pi b ps];
var_select = linspace(1,M,M);
nH = length(vd_interval);
hendval = 1;
hrz = [linspace(0,hendval,nH-1)'; hendval*1.2];
Cell{1} = {num2str(vd_interval(1:end-1)-1);'\infty'};
y_nm_full = ["output";"consumption";"investment";"wage";"labor";"market value of debt";"primary surplus";"govt spending";"transfer";"inflation";"interest rate";"nominal bond return";"output gap";"primary surplus"];
y_nm = y_nm_full(data_index);
mycol = {'#00539a' '#ff832b' '#198038' '#a56eff'};
fsz_tl = 18;
ftsz2 = 18;
%%
figure('name','Variance Decomposition: MP shock','color','w','WindowState','maximized')
if length(var_select) == 9
    tiledlayout(3,3, 'Padding', 'none', 'TileSpacing', 'compact');    
elseif length(var_select) == 8
    tiledlayout(2,length(var_select)/2, 'Padding', 'none', 'TileSpacing', 'compact');    
else
    tiledlayout(1,length(var_select), 'Padding', 'none', 'TileSpacing', 'compact');    
end
for j = var_select
    nexttile
    bar(hrz,[mean(FEVD(j,vd_interval,:),3)' mean(FEVD_DSGE(j,vd_interval,:),3)'],1,'edgecolor','w')    
    colororder(mycol)
    %axis tight    
    ticklabel = get(gca,'TickLabel');
    set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz2,'FontWeight','bold');
    set(gca, 'XTick',hrz, 'XTickLabel',Cell{1})    
    title(y_nm(j),'FontName','Consolas','fontsize',fsz_tl,'FontWeight','bold');
    if j == 1
       legend('VAR','DSGE','location','northwest','FontSize',ftsz2,'FontName','Consolas','FontWeight','bold');
       legend boxoff 
    end    
end