function [IRF DVD FEVD] = sign_restriction(BVAR_est,H,sign_irf,shock_index)

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

% measure the sign restriction horizon
[~, H_sign] = size(sign_irf);

% IRF and Debt value decomposition basket
IRF = zeros(M,H+1,N);
FETV = zeros(M,H+1,N);
IDSV = zeros(M,H+1,N);
DVD = zeros(N,4);
Q_draw_trial = 0;

% time discount rate for Debt value decomposition
betat = .99;

% basis vectors for variable selection
cps = double(1:M*p == ps)';
cpi = double(1:M*p == Pi)';
cRB = double(1:M*p == RB)';
cb = double(1:M*p == b)';

% draw IRF
for j = 1:N
    % retrieve coefficient matrix from vectorization
    Aj = A_sampler(:,:,j);
    
    % is there a constant vector in VAR?
    if regressor_constant
        Aj = Aj(2:end,:)';
    else
        Aj = Aj';
    end
    
    % companion VAR
    Aj_comp = [Aj; [kron(eye(p-1),eye(M)) zeros(M*(p-1),M)]];
    
    % draw Sigma from sampler
    Sigmaj = Sigma_sampler(:,:,j);
    
    % cholesky decomposition
    Lj = chol(Sigmaj,'lower');
    
    % starting ticket of sign restriction
    sign_uniqueness = 1e2;
    sign_match = 0;
    
    % sign restriction process
    while sign_uniqueness > 1 || sign_match < 1
        
        % count up the Q draw trial
        Q_draw_trial = Q_draw_trial + 1;
        
        % QR decomposition to draw sign matrix from Haars distribution
        QR = randn(M,M);
        [Q,~] = qr(QR);
        
        % check the uniqueness of the sign vector
        % Wolf(2019) points out the problem here...
        Signj = sign(Lj*Q);
        sign_uniqueness = sum(prod(Signj(:,shock_index) == Signj));
        
        % compute IRF upto the sign restricted horizon
        irfj = [];
        fetv = Sigmaj;
        Bj = Lj*Q*double(1:M==shock_index)';
        Bj = Bj / Bj(R) * .25;
        for h = 0:H_sign-1
            Ah = Aj_comp^h;
            irfj = [irfj Ah(1:M,1:M)*Bj];
            FETV(:,h+1,j) = diag(Ah(1:M,1:M)*Sigmaj*Ah(1:M,1:M)');
            IDSV(:,h+1,j) = diag(Ah(1:M,1:M)*Bj*Bj'*Ah(1:M,1:M)');
        end
        
        % evaluate the sign
        if sum(sum(sign_irf .* sign(irfj) < 0)) == 0
            IRF(:,1:H_sign,j) = irfj;
            impact = Lj*Q*double(1:M==shock_index)';            
            sign_match = sign_match + 1;
            % perhaps the opposite may work!
        elseif sum(sum(sign_irf .* sign(-irfj) < 0)) == 0
            IRF(:,1:H_sign,j) = -irfj;
            impact = -Lj*Q*double(1:M==shock_index)';            
            sign_match = sign_match + 1;
        end        
    end    
    
    % compute the IRF beyond the restricted horizon
    irfj = [];
    impact = impact / impact(R) * .25;
    for h = H_sign:H
        Ah = Aj_comp^h;
        irfj = [irfj Ah(1:M,1:M)*impact];
        FETV(:,h+1,j) = diag(Ah(1:M,1:M)*Sigmaj*Ah(1:M,1:M)');
        IDSV(:,h+1,j) = diag(Ah(1:M,1:M)*impact*impact'*Ah(1:M,1:M)');
    end
    IRF(:,H_sign+1:end,j) = irfj;
    
    
    % impact vector for debt value decomposition
    impact_vec = [impact; zeros(M*(p-1),1)];
    
    % debt value decomposition
    inflation_impact = cpi' * impact_vec;
    bond_price_reval = cRB' * impact_vec;
    debt_reval = bond_price_reval - inflation_impact;
    Ex_dc_path = -(cRB - cpi)' * inv(eye(M*p) - betat * Aj_comp) * impact_vec + debt_reval;
    Ex_ps_path = cps' * inv(eye(M*p) - betat * Aj_comp) * impact_vec;
    
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