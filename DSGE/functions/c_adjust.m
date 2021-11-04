function c = c_adjust(x,z)
% scaling factor adjustment with target acceptance rate z
% x is the last acceptance rate
% if one wants to use constant factor, set x = z;
A = .95;
B = 2*(1-A);
D = 1;
c = A + B * exp(D*(x-z)) / (1 + exp(D*(x-z)));
end