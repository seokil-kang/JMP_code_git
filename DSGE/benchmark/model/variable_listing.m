endogenous_var_list = {
    % allocations
    'c' % cousumption
    'ivst' % investment
    'n' % labor
    'v' % capital utilization
    'b' % long-term govt bond
    'y' % output
    'k' % utilized capital
    'ks' % capital stock
    
    % prices
    'infl' % (pi) inflation
    'w' % real wage
    'rk' % capital rental rate
    'R' % gross nominal interest rate
    'Q' % nominal long-term govt bond price
    'RL' % nominal long-term govt bond return
    
    % govt policy
    'g' % govt spending
    'z' % lump-sum transfer
    'tau' % total tax revenue
    's' % primary surplus
    
    % proxy variables
    'Uc' % marginal utility of consumption
    'lambda' % Lagrange multiplier on HHBC
    'wmk' % wage markup
    'pmk' % wage markup
    'q' % Tobin's q
    'mktb' % market value of debt
    };
    
    % exogenous variables
    shock_var_list = {
    'u_gamma' % (gamma) trend growth(labor-augmented permanent technology growth)
    'u_U' % int-temp preference
    'u_i' % (normalized) marginal efficiency of investment
    'u_w' % (normalized) wage markup
    'u_p' % (normalized) price markup
    'u_g' % govt spending shock
    'u_z' % lump-sum transfer shock
    'u_M' % MP shock
    };

    obs_var_list = {
    };

    lag_var_list = {
    'Ly'
    'Livst'
    'Lw'
    'Lb'
    };

    expectation_var_list = {
    % expectation terms
    'xu_gamma'
    'xUc'
    'xiv'
    'xlambda'
    'xrk'
    'xq'
    'xRL'
    'xinfl'
    'xw'
    };

    exogenous_var_list = {
    };

    flexible_var_list = {
    };

    var_name_list = [
                    endogenous_var_list;
                    obs_var_list;
                    lag_var_list;
                    expectation_var_list;                    
                    flexible_var_list
                    shock_var_list
                    ];
    n_var = length(var_name_list);
    n_shock = length(shock_var_list);
    n_forecast_error = length(expectation_var_list);
    for var_index = 1:n_var
        eval([char(var_name_list(var_index)) ' = var_index;']);
        eval([char(strcat({'var_name.'},var_name_list(var_index))) ' = var_index;']);
    end
    for forecast_error_index = 1:n_forecast_error
        eval([char(strcat({'eta_'},expectation_var_list(forecast_error_index))) ' = forecast_error_index;']);
    end
    for shock_index = 1:n_shock
        eval([char(strcat({'eps_'},shock_var_list(shock_index))) ' = shock_index;']);
    end