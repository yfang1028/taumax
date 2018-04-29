function [x, tau, c, m, w, iter] = optimizex6_acor(uall, CnM, x0, M)
warning ("off");

x = x0;
xprev = 0;
lam = 1;
lamprev = 0;
N = size(uall, 2);
iter = 0;
c = zeros(1, M+1);

m = 1;

while norm(x-xprev)>1e-6 && iter<1e2
  u = x'*uall;

  [c, uncf] = acffft(u, N);
  lamprev = lam;
  w = lagwindowacor(c', M+1);
  m = sum(w);

  xprev = x;

  T = T + sum(permute(w(1:M).*permute(CnM(:,:,2:M+1),[1,3,2]),[1,3,2]),3);

  [v, eigvs] = eig(T,CnM(:,:,1)); 
  [cc, idx] = max(abs(diag(eigvs)));
  x = v(:, idx);
  iter = iter+1;
  
end

tau = sum(c(2:m))*2+1;



end