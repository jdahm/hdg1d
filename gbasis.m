function GPhi = gbasis(xn, x);
% evaluates basis function gradients at x

n = length(x);
order = length(xn)-1;
GPhi = zeros(n, order+1);
for p=0:order,
  B = glagrange(xn, p+1,x);
  GPhi(:,p+1) = B';
end