function integration_test()

% f = x^p
p = 6
[x,w] = lgwt(ceil((p+1)/2), 0, 1);
y = x.^p;
sum(w.*y) - 1.0/(p+1)

