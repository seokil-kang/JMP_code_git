function IR = impulse_response(G,M,H)
% impulse response of the state space model
% X(t) = C + G * X(t-1) + M * e_mp;
[n_var n_shock] = size(M);
IR = zeros(n_var,H+1,n_shock);
for h = 0:H
    IR(:,h+1,:) = G^h * M;
end

end