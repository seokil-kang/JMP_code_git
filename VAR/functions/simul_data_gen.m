function [y A_comp Ac IRF L] = simul_data_gen(N,M,H)
burnin = 5e3;
p = 4;
initial_var = randn(M,4);
A1 = eye(M)*.2 + diag(rand(M,M))*.02 + rand(M,M)*0.05;
A2 = eye(M)*.15 + diag(rand(M,M))*.01 + randn(M,M)*0.05;
A3 = eye(M)*.05 + diag(rand(M,M))*.01 + randn(M,M)*0.01;
A4 = eye(M)*.01 + diag(rand(M,M))*.005 + randn(M,M)*0.01;
Ac = randn(M,1)*0.01;
%Ac(randperm(M,3)) = 2;
L = randn(M,M);
[L ~] = qr(L);
y1 = Ac + A1*initial_var(:,4) + A2*initial_var(:,3) + A3*initial_var(:,2) + A4*initial_var(:,1) + L*randn(M,1);
y2 = Ac + A1*y1 + A2*initial_var(:,4) + A3*initial_var(:,3) + A4*initial_var(:,2) + L*randn(M,1);
y3 = Ac + A1*y2 + A2*y1 + A3*initial_var(:,4) + A4*initial_var(:,3) + L*randn(M,1);
y4 = Ac + A1*y3 + A2*y2 + A3*y1 + A4*initial_var(:,4) + L*randn(M,1);
y5 = Ac + A1*y4 + A2*y3 + A3*y2 + A4*y1 + L*randn(M,1);
y = [y1 y2 y3 y4 y5];
for j = 1:N+burnin
    yt = Ac + A1*y(:,end) + A2*y(:,end-1) + A3*y(:,end-2) + A4*y(:,end-3) + L*randn(M,1);
    y = [y yt];
end
y = y(:,end-N+1:end)';

A_comp = [A1 A2 A3 A4; [kron(eye(p-1),eye(M)) zeros(M*(p-1),M)]];

IRF = [];
impact = zeros(M,1);
impact(1) = 1;
impact = L*impact;
for h = 0:H
    Ah = A_comp^h;
    IRF = [IRF Ah(1:M,1:M)*impact];
end
IRF = IRF';
hrz = linspace(0,H,H+1)';
if 0
    for j = 1:M
        subplot(3,3,j)
        hold on
        yline(0,'r')
        plot(hrz,IRF(:,j),'x-','linewidth',1.5)
        hold off
        xlim([0 H])
    end
end