function [X_smooth X_update G C M] = kalman_smoother(theta,Y,ME_scale)
%two-sided Kalman filter for the following linear gaussian state space model
%transition eqn: X(t) = C + G * X(t-1) + M * e(t), e ~ N(0,S)
%measurement eqn: Y(t) = D + H * X(t) + W * u(t), u ~ N(0,V)
%where V = ME_scale*var(Y)
% - inputs
%(1) theta: model parameter.
%(2) Y: observations
%(3) ME_scale: measurement error scale. If 0, no measurement error

% - output
%X_KS: Kalman Smoothed state vector

%- note: Using mvnpdf is not recommended, as it requires Positive Semi Definite matrix.
%        Numerical process usually does not preserve it. Correcting it slows down the process.

% computes boundary conditions first
boundcheck = preliminary_condition(theta);

% there is no need to compute likelihood if the parameter fails boundary conditions
if boundcheck == 0
    X_smooth = nan;
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
        V = ME_scale*V;
        
        % measure the state space model dimension & observables
        n_var = length(C);
        [T, n_obs] = size(Y);
        
        % kalman filter initializing (L means lag operator)
        
        % 1st moment for state vector: long-run mean at t = 0
        Xtt = inv(eye(n_var)-G)*C;
        
        % 2nd moment from Lyapnouv equation at t = 0
        SIGtt = dlyap(G, M * S * M');
        
        % filter basket
        X_predict = [];
        X_update = [Xtt];
        SIG_predict = [];
        SIG_update = [SIGtt];        
        
        for t = 1:T
            
            % yesterday's t is today's Lt
            XLtLt = Xtt;
            SIGLtLt = SIGtt;
            
            %(a) prediction step
            
            % 1st moment prediction
            XtLt = C + G * XLtLt;
            YtLt = D + H * XtLt;
            
            % 2nd moment prediction
            SIGtLt = M * S * M' + G * SIGLtLt * G';
            Omega_t = W * V * W' + H * SIGtLt * H';
            
            % collect the prediction
            X_predict = [X_predict XtLt];
            SIG_predict = [SIG_predict SIGtLt];
            
            % determinant for covariance matrix Omega
            detOmega_t = det(Omega_t);
            
            % its determinant must be weekly positive
            if detOmega_t < 0
                X_smooth = nan;
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
            epsilon_t = Y_t - YtLt;
                        
            %(c) Updating step
            
            % compute Kalman gain
            K_t = SIGtLt * H' * invOmega_t;
            
            % update 1st moment
            Xtt = XtLt + K_t * epsilon_t;
            
            % update 2nd moment
            SIGtt = SIGtLt - K_t * H * SIGtLt;
            
            % collect the update
            X_update = [X_update Xtt];
            SIG_update = [SIG_update SIGtt];
        end
        
        % compute Kalman smoother backward
        X_smooth = [Xtt];
        % use pseudo-inverse
        % Too small tolerance gives Jt = 0
        % Too large tolerance gives X_smooth = X_update
        for tb = 1:T
            tol = 1e-9;
            Jt = SIG_update(:,end-(tb+1)*n_var+1:end-tb*n_var)*G'*pinv(SIG_predict(:,end-tb*n_var+1:end-(tb-1)*n_var),tol);
            X_smooth = [X_update(:,end-tb) + Jt*(X_smooth(:,end)-X_predict(:,end-tb+1)) X_smooth];
        end            
        
    else % if no unique&stable solution
        X_smooth = nan;
    end
end
end