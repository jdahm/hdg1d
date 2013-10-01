function gphi = glagrange(xn, j, x);
% computes values of the gradient of the jth lagrange function 
% based on nodes in xn, at the x-values given in x

n = length(xn);
nj = [1:j-1,j+1:n]';
den = prod(xn(j)-xn(nj));
if (length(nj) == 0)
  gphi = zeros(size(x)); return;
elseif (length(nj) == 1)
  num = ones(length(x),1);
else
  num = zeros(length(x),1);
  for i=nj',
    nij = setdiff(nj,i);
    xnij = xn(nij);
    num = num + prod(repmat(x,1,n-2) - repmat(xnij',length(x),1), 2);
  end
end
gphi = num/den;