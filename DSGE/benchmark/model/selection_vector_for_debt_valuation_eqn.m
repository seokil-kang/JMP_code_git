% construct selection vectors for debt valuation response
function [cRL cINFL cPS cMP] = selection_vector_for_debt_valuation_eqn

% call variable indexing
variable_listing

% selection vector for nominal return of govt bond
cRL = zeros(n_var,1);
cRL(RL) = 1;

% selection vector for inflation
cINFL = zeros(n_var,1);
cINFL(infl) = 1;

% selection vector for primary surplus
cPS = zeros(n_var,1);
cPS(s) = 1;

% selection vector for monetary policy shock
cMP = zeros(n_shock,1);
cMP(eps_u_M) = 1;

end