function F = force_mix_gauss_2d(x)

d = size(x,1);

lam1 = [36;1]; lam2 = [1;36]; mu1 = [-1;3]; mu2 = [2;0];
rho1 = exp(-1/2*(x-mu1)'*(lam1.*(x-mu1)));
rho2 = exp(-1/2*(x-mu2)'*(lam2.*(x-mu2)));
Ux = (lam1.*(x-mu1)*rho1 + lam2.*(x-mu2)*rho2)/(rho1 + rho2);

F = -Ux;

end