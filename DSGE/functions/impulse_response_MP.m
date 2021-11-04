function IR = impulse_response_MP(G,M,H)
% impulse response of the state space model
% X(t) = C + G * X(t-1) + M * e_mp;
[n_var ~] = size(M);
IR = zeros(n_var,H+1);
for h = 0:H
    IR(:,h+1) = G^h * M(:,end);
end

end