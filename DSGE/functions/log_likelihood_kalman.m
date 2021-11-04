function ln_likelihood = log_likelihood_kalman(theta,Y,ME_scale)
%Compute model likelihood with Kalman filtering
%for the following linear gaussian state space model
%transition eqn: X(t) = C + G * X(t-1) + M * e(t), e ~ N(0,S)
%measurement eqn: Y(t) = D + H * X(t) + W * u(t), u ~ N(0,V)
%where V = ME_scale*var(Y)
% - inputs
%(1) theta: model parameter.
%(2) Y: observations
%(3) ME_scale: measurement error scale. If 0, no measurement error

% - output
%ln_likelihood: log likelihood of model.

%- note: Using mvnpdf is not recommended, as it requires Positive Semi Definite matrix.
%        Numerical process usually does not preserve it. Correcting it slows down the process.

% computes boundary conditions first
boundcheck = preliminary_condition(theta);

% there is no need to compute likelihood if the parameter fails boundary conditions
if boundcheck == 0
    ln_likelihood = -inf;
else % proceed to filtering
    
    % Get the model equilibrium condition system
    [G0 G1 CC Psi Pi S] = equilibrium_condition(theta);
    
    % solve the model
    [G,C,M,~,~,~,~,eu] = gensys(G0,G1,CC,Psi,Pi);
    
    % evaluate likelihood for unique and stable solutions only
    if (eu(1) == 1) && (eu(2) == 1)
        
        % get mesurement equation system
        [D H W V] = measurement_equation(theta);
        
        % measurement error specification
        %V = ME_scale*cov(Y);
        V = ME_scale*V;
        
        % measure the state space model dimension & observables
        n_var = length(C);
        [T, n_obs] = size(Y);
        
        % kalman filter initializing (L means lag operator)
        
        % 1st moment for state vector: long-run mean at t = 0
        EtXt = inv(eye(n_var)-G)*C;
        
        % 2nd moment from Lyapnouv equation at t = 0
        EtSIGt = dlyap(G, M * S * M');
        
        % initial log likelihood
        ln_likelihood = 0;
        
        for t = 1:T
            
            % yesterday's t is today's Lt
            ELtXLt = EtXt;
            ELtSIGLt = EtSIGt;
            
            %(a) prediction step
            
            % 1st moment prediction
            ELtXt = C + G * ELtXLt;
            ELtYt = D + H * ELtXt;
            
            % 2nd moment prediction
            ELtSIGt = M * S * M' + G * ELtSIGLt * G';
            Omega_t = W * V * W' + H * ELtSIGt * H';
            
            % determinant for covariance matrix Omega
            detOmega_t = det(Omega_t);
            
            % its determinant must be weekly positive
            if detOmega_t < 0
                ln_likelihood = -inf;
                display('Filtered covariance has negative determinant')
                return
            end
            
            % compute the inverse of Omega
            invOmega_t = inv(Omega_t);
            
            %(b) likelihood evaluation
            
            % pick up the observable set at t and transpose it to column
            Y_t = Y(t,:)';
            
            % filtered likelihood follows normal distriubtion with
            % mean ELtYt and covariance term Omega_t
            
            % residual of filtering
            epsilon_t = Y_t - ELtYt;
            
            % likelihood at t
            lnl_t = -0.5 * n_obs * log(2*pi) - 0.5 * log(detOmega_t)...
                -0.5 * epsilon_t' * invOmega_t * epsilon_t;
            
            % add up the log likelihood function repeatedly
            ln_likelihood = ln_likelihood + lnl_t;
            
            %(c) Updating step
            
            % compute Kalman gain
            K_t = ELtSIGt * H' * invOmega_t;
            
            % update 1st moment
            EtXt = ELtXt + K_t * epsilon_t;
            
            % update 2nd moment
            EtSIGt = ELtSIGt - K_t * H * ELtSIGt;
        end
        
    else % if no unique&stable solution
        ln_likelihood = -inf;
    end
end
end