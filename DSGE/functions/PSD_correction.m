function Y = PSD_correction(X)
if issymmetric(X)&& min(eig(X))> 0
    Y = X;
else
if issymmetric(X)
else
    X = (X+X')/2;
end
[Lchol, DMC, Ptb] = modchol_ldlt(X);
PL = Ptb'*Lchol*DMC^(1/2);
Y = PL*PL';
end
end