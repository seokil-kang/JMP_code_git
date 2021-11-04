function [IRF DVD FEVD Phi] = proxy_iv(Y,X,BVAR_est,H,Z,policy_index)

% unwrap the structure
A_sampler = BVAR_est.A_sampler;
Sigma_sampler = BVAR_est.Sigma_sampler;
data_index = BVAR_est.data_index;
diff_index = BVAR_est.diff_index;

% call the variable indexing
variable_indexing(data_index);

% measure the VAR dimension size
[M,~,~] = size(Sigma_sampler);

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

% IRF and Debt value decomposition basket
IRF = zeros(M,H+1,N);
FETV = zeros(M,H+1,N);
IDSV = zeros(M,H+1,N);
DVD = zeros(N,4);
Phi = [];

% time discount rate for Debt value decomposition
betat = .99;

% basis vectors for variable selection
cps = double(1:M*p == ps)';
cpi = double(1:M*p == Pi)';
cRB = double(1:M*p == RB)';
cb = double(1:M*p == b)';

% non policy indice
q_index = setdiff(1:M,policy_index);

% baskets to discard
j_no_count = [];

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
    
    %{
    % subblocks of the covariance matrix
    Sigma11 = Sigmaj_ro(policy_index,policy_index);
    Sigma12 = Sigmaj_ro(policy_index,q_index);
    Sigma21 = Sigma12';
    Sigma22 = Sigmaj_ro(q_index,q_index);
    
    % IV estimation (notation follows Lakdawala(JAE, 2019))
    % step 1: estimate B21B11^-1
    B21invB11 = (inv(Z' * upj) * (Z' * uqj))';
    
    % step 2: some large proxy matrix for simplicity
    Q = B21invB11 * Sigma11 * B21invB11' - (Sigma21 * B21invB11' + B21invB11 * Sigma12) + Sigma22;
    
    % step 3: compute B12B12'
    B12B12T = (Sigma21 - B21invB11 * Sigma11)' * inv(Q) * (Sigma21 - B21invB11 * Sigma11);
    
    % step 4: compute B11B11'
    B11B11T = Sigma11 - B12B12T;
    
    % discard negative draws
    if B11B11T < 0
        j_no_count = [j_no_count j];
    else        
        % scalar makes life easier
        B11 = sqrt(B11B11T);
        
        % step 5: compute B21
        B21 = B21invB11 * B11;
        
        % identified impact of the policy shock
        impact = zeros(M,1);
        impact(policy_index) = B11;
        impact(q_index) = B21;
        
        % normalize to 25 basis-points shock
        %impact = impact / impact(policy_index) * .25;
        
        % is there a constant vector in VAR?
        if regressor_constant
            Aj = Aj(2:end,:)';
        else
            Aj = Aj';
        end
        
        % companion VAR
        Aj_comp = [Aj; [kron(eye(p-1),eye(M)) zeros(M*(p-1),M)]];
        
        % compute the IRF beyond the restricted horizon
        irfj = [];
        for h = 0:H
            Ah = Aj_comp^h;
            irfj = [irfj Ah(1:M,1:M)*impact];
        end
        IRF(:,:,j) = irfj;
        
        % impact vector for debt value decomposition
        impact_vec = [impact; zeros(M*(p-1),1)];
        
        % debt value decomposition
        inflation_impact = cpi' * impact_vec;
        bond_price_reval = cRB' * impact_vec;
        debt_reval = bond_price_reval - inflation_impact;
        Ex_dc_path = -(cRB - cpi)' * inv(eye(M*p) - betat * Aj_comp) * impact_vec + debt_reval;
        Ex_ps_path = (1 - betat) * cps' * inv(eye(M*p) - betat * Aj_comp) * impact_vec;
        
        % collect the decomposition
        DVD(j,:) = [inflation_impact bond_price_reval Ex_dc_path Ex_ps_path];  
    end
    %}
    
    % A nice approach for univariate proxy IV case(Lunsford 2015)
    
    % measure the sample size for IV estimation
    T = length(Z);
    
    % recover the IV correlation
    phij = sqrt((Z'*Uj/T)*inv(Uj'*Uj)*(Uj'*Z));
    Phi = [Phi; phij];
    
    % estimate the impact
    impact = (Uj'*Z/T)/phij;
    
    % contractionary shock
    if impact(policy_index) < 0
        impact = -impact;
    end
    
    % normalize to 25 basis-points shock
    %impact = impact / impact(policy_index) * .25;
    
    % is there a constant vector in VAR?
    if regressor_constant
        Aj = Aj(2:end,:)';
    else
        Aj = Aj';
    end
    
    % companion VAR
    Aj_comp = [Aj; [kron(eye(p-1),eye(M)) zeros(M*(p-1),M)]];
    
    % companion impact vector
    impact_comp = [impact; zeros(M*(p-1),1)];
    
    % compute the IRF beyond the restricted horizon
    irfj = impact_comp;
    forcast_error_total_variance = kron(reshape(double(1:p*p == 1),p,p),Sigmaj);
    ID_shock_variance = impact_comp * impact_comp';
    for h = 0:H
        IRF(:,h+1,j) = irfj(1:M,:);
        irfj = Aj_comp*irfj;
        FETV(:,h+1,j) = diag(forcast_error_total_variance(1:M,1:M));
        IDSV(:,h+1,j) = diag(ID_shock_variance(1:M,1:M));
        forcast_error_total_variance = Aj_comp*forcast_error_total_variance*Aj_comp';
        ID_shock_variance = Aj_comp * ID_shock_variance * Aj_comp';
    end    
        
    % debt value decomposition
    inflation_impact = cpi' * impact_comp;
    bond_price_reval = cRB' * impact_comp;
    debt_reval = bond_price_reval - inflation_impact;
    Ex_dc_path = -(cRB - cpi)' * inv(eye(M*p) - betat * Aj_comp) * impact_comp + debt_reval;
    Ex_ps_path = cps' * inv(eye(M*p) - betat * Aj_comp) * impact_comp;
    
    % collect the decomposition
    DVD(j,:) = [inflation_impact bond_price_reval Ex_dc_path Ex_ps_path];  
    
end

% cumulate the IRF for differenced data
if length(diff_index) > 0
    IRF(diff_index,:,:) = cumsum(IRF(diff_index,:,:),2);
end

% Forecast Error decomposition
FEVD = cumsum(IDSV,2) ./ cumsum(FETV,2);
end