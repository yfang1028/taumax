function m = optm(lam, c0, sigma2)

sigma2 = sigma2/c0^2;

a = 2*lam^2/(1-lam^2)^2+2*sigma2;
b = 2*sigma2*(1+lam);
c = sigma2/log(lam);

mu = (b+sqrt(b^2-4*a*c))/2/a;
m = log(mu)/log(lam);
m = floor(m);