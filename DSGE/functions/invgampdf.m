function y = invgampdf(x,a,b)
% inverse gamma probability density function based on gamma pdf.
% X ~ G(a,1/b) <=> Y = 1/X ~ IG(a,b)
% PY(y) = PX(x) / y^2
y = gampdf(1./x,a,1/b)./(x.^2);
end