% compute the IRF & debt valuation decompositon to contractionary MP shock

function [IRF DVD FEVD] = IRF_and_DVD(theta,H)

% add path for submodules
%addpath('functions','gensys','model','prior','sims_minwel')

% measure some sizes
[n,N] = size(theta);

% call the variable labels
variable_listing;

% selection vectors
cRL = double(1:n_var == RL)';
cinfl = double(1:n_var == infl)';
cs = double(1:n_var == s)';
cMP = double(1:n_shock == eps_u_M)';

% baskets for computation results
IRF = [];
DVD = [];
FEVD = [];

% main computation
for j = 1:N
    
    % select parameter set
    thetaj = theta(:,j);
    
    % setup the model
    [G0 G1 C Psi Pi] = equilibrium_condition(thetaj);
    
    % solve the model
    [G,C,M,~,~,~,~,eu] = gensys(G0,G1,C,Psi,Pi);
    
    % normalize to 25 bps.
    %M = M/M(R,eps_u_M)*.25;
    
    % compute only if we have uniqueness and stability
    if (eu(1) == 1) && (eu(2) == 1)
        
        % IRF computation
        irfj = [];
        fetv = [];
        idsv = [];
        for h = 0:H            
            irfj = [irfj G^h * M(:,eps_u_M)];
            if h <= min(H,40)
                fetv = [fetv diag(G^h*M*(G^h*M)')];
            end
        end
        for m = 1:n_shock
            idsvm = [];
            for h = 0:min(H,40)
                idsvm = [idsvm diag(G^h*M(:,m)*(G^h*M(:,m))')];
            end
            idsv = cat(3,idsv,cumsum(idsvm,2)./cumsum(fetv,2));
        end
        
        
        % stackup the results in 3D
        IRF = cat(3,IRF,irfj);
        FEVD = cat(4,FEVD,idsv);
        
        % compute some steady states and functions of parameters
        [betat bss] = debt_valuation_response_setup(thetaj);
        
        % Debt valuation factors
        inflation_impact = cinfl' * M * cMP;
        bond_price_reval = cRL' * M * cMP;
        debt_reval = bond_price_reval - inflation_impact;
        Ex_dc_path = -(cRL - cinfl)' * inv(eye(n_var) - betat * G) * M * cMP + debt_reval;
        Ex_ps_path = cs' * inv(eye(n_var) - betat * G) * M * cMP;
        
        % stackup the results
        DVD = [DVD; inflation_impact bond_price_reval Ex_dc_path Ex_ps_path];
    end
end
end