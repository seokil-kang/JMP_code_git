function [X nm nmtx] = prior_info_setup
% setup the prior information of each parameter
% 1st column: distribution calirfication
% 1: Beta, 2: Gamma, 3: Inverse-Gamma, 4: Normal, 5: Uniform
% 2nd column: Mean of prior distribution
% 3rd column: Std dev of prior distribution
X = [];
nm = [];
nmtx = [];

% permanent growth rate *100
X = [X; 4, 0.4, 0.05];
nm = [nm; "gamm100"];
nmtx = [nmtx; "$100\gamma$"];

% constant relative degree of risk aversion
%X = [X; 2, 1.5, 0.5];
%nm = [nm; "chi_c"];
%nmtx = [nmtx; "$\chi_n$"];

% inverse Frisch labor elasticity
X = [X; 2, 2, 0.5];
nm = [nm; "chi_n"];
nmtx = [nmtx; "$\chi_n$"];

% habit formation
X = [X; 1, 0.5, 0.2];
nm = [nm; "eta"];
nmtx = [nmtx; "$\eta$"];

% capital utilization parameter
X = [X; 1, 0.6, 0.15];
nm = [nm; "nu"];
nmtx = [nmtx; "$\nu$"];

% marginal cost of investment adjustment
X = [X; 4, 6, 1.5];
nm = [nm; "varsigma"];
nmtx = [nmtx; "$\varsigma$"];

% price rigidity
X = [X; 1, 0.5, 0.2];
nm = [nm; "zeta_p"];
nmtx = [nmtx; "$\zeta_p$"];

% wage rigidity
X = [X; 1, 0.5, 0.2];
nm = [nm; "zeta_w"];
nmtx = [nmtx; "$\zeta_w$"];

% price indexation
X = [X; 1, 0.5, 0.2];
nm = [nm; "iota_p"];
nmtx = [nmtx; "$\iota_p$"];

% wage indexation
X = [X; 1, 0.5, 0.2];
nm = [nm; "iota_w"];
nmtx = [nmtx; "$\iota_w$"];

% MP persistency
X = [X; 1, 0.5, 0.2];
nm = [nm; "rho_R"];
nmtx = [nmtx; "$\rho_R$"];

% MP response to inflation
X = [X; 5, 2.04, 0.6];
nm = [nm; "phi_pi"];
nmtx = [nmtx; "$\phi_\pi$"];

% MP response to output gap
X = [X; 4, 0.15, 0.2];
nm = [nm; "phi_x"];
nmtx = [nmtx; "$\phi_x$"];

% FP(gov't spending) persistency
X = [X; 1, .5, .2];
nm = [nm; "rho_g"];
nmtx = [nmtx; "$\rho_g$"];

% FP(gov't spending) response to output gap
X = [X; 4, 0.15, 0.5];
nm = [nm; "varphi_x"];
nmtx = [nmtx; "$\varphi_x$"];

% FP(gov't spending) response to debt
X = [X; 4, 1.2, 0.5];
nm = [nm; "varphi_b"];
nmtx = [nmtx; "$\varphi_b$"];

% FP(transfer) response to output gap
X = [X; 4, 0.15, 0.5];
nm = [nm; "psi_x"];
nmtx = [nmtx; "$\psi_x$"];

% FP(transfer) response to debt
X = [X; 4, 1.2, 0.5];
nm = [nm; "psi_b"];
nmtx = [nmtx; "$\psi_b$"];

% permanent growth persistency
X = [X; 1, 0.5, 0.2];
nm = [nm; "varrho_gamma"];
nmtx = [nmtx; "$\varrho_\gamma$"];

% preference shock persistency
X = [X; 1, 0.5, 0.2];
nm = [nm; "varrho_U"];
nmtx = [nmtx; "$\varrho_U$"];

% marginal investment efficiency persistency
X = [X; 1, 0.5, 0.2];
nm = [nm; "varrho_i"];
nmtx = [nmtx; "$\varrho_i$"];

% wage markup persistency
X = [X; 1, 0.5, 0.2];
nm = [nm; "varrho_w"];
nmtx = [nmtx; "$\varrho_w$"];

% price markup persistency
X = [X; 1, 0.5, 0.2];
nm = [nm; "varrho_p"];
nmtx = [nmtx; "$\varrho_p$"];

% gov't spending persistency
X = [X; 1, 0.5, 0.2];
nm = [nm; "varrho_g"];
nmtx = [nmtx; "$\varrho_g$"];

% MP shock persistency
X = [X; 1, 0.5, 0.2];
nm = [nm; "varrho_M"];
nmtx = [nmtx; "$\varrho_M$"];

% permanent growth std dev *100
X = [X; 3, 0.1, 1];
nm = [nm; "sigma_gamma"];
nmtx = [nmtx; "$\sigma_\gamma$"];

% preference shock std dev *100
X = [X; 3, 0.1, 1];
nm = [nm; "sigma_U"];
nmtx = [nmtx; "$\sigma_U$"];

% marginal investment efficiency std dev *100
X = [X; 3, 0.1, 1];
nm = [nm; "sigma_i"];
nmtx = [nmtx; "$\sigma_i$"];

% wage markup std dev *100
X = [X; 3, 0.1, 1];
nm = [nm; "sigma_w"];
nmtx = [nmtx; "$\sigma_w$"];

% price markup std dev *100
X = [X; 3, 0.1, 1];
nm = [nm; "sigma_p"];
nmtx = [nmtx; "$\sigma_p$"];

% gov't spending std dev *100
X = [X; 3, 0.1, 1];
nm = [nm; "sigma_g"];
nmtx = [nmtx; "$\sigma_g$"];

% transfer shock std dev *100
X = [X; 3, 0.1, 1];
nm = [nm; "sigma_z"];
nmtx = [nmtx; "$\sigma_z$"];

% MP shock std dev *100
X = [X; 3, 0.1, 1];
nm = [nm; "sigma_M"];
nmtx = [nmtx; "$\sigma_M$"];

% mean of measured hours work
X = [X; 4, 0, 1];
nm = [nm; "n_bar"];
nmtx = [nmtx; "$\bar{n}$"];

% mean of measured inflation
X = [X; 4, 0.79, 0.25];
nm = [nm; "pi_bar"];
nmtx = [nmtx; "$\bar{\pi}$"];

% mean of measured primary surplus
X = [X; 4, .21, .5];
nm = [nm; "s_bar"];
nmtx = [nmtx; "$\bar{s}$"];

% std dev of measurement error for primary surplus
X = [X; 3, .1, 1];
nm = [nm; "sigma_ME_s"];
nmtx = [nmtx; "$\sigma_{ME_s}$"];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save the parameter name to the model folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('model/parameter_name','nm');

% prior distribution specification
prior_spec = X(:,1);
save('model/prior_spec','prior_spec');

end