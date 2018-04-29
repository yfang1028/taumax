function [x, tau, c, u, m, w, lam, c0, v, iter, isfail] = optimizex6(uall, CnM, x0, M)
warning ("off");
format long;
MINESS = 100;

isfail = 0;
x = x0;
xprev = 0;
lam = 1;
lamprev = 0;
N = size(uall, 2);

iter = 0;
c = zeros(1, M+1);
u = zeros(1, M+1);
c0 = 1;
m(1) = 1;

while norm(lam-lamprev)>1e-6 && iter<1e2

  u = x'*uall;
  cprev = c;
  c0prev = c0;
  uprev = u;
  [c, uncf] = acffft(u, N);
  lamprev = lam;
  [lam, c0, sigma2, ii] = findlam2(uncf, 0, M); 
  c = uncf;
  
  if lam > 1
    c = cprev;
    c0 = c0prev;
    lam = lamprev;
    u = uprev;
  end
  if iter > 0 && lam < lamprev
    c = cprev;
    c0 = c0prev;
    lam = lamprev;
    u = uprev;
  end
    
  if (1+lam)/(1-lam) > N/MINESS
    isfail = 1;
  end
  

  m1 = optm(lam, c0, sigma2);
  m(iter+2) = m1;
  w = lagwindowours(m1, M-1, lam);
  w = w';
  if length(w) > M-1
    w = w(1:M-1);
  end

  xprev = x;

  T = T + sum(permute(w'.*permute(CnM(:,:,2:M),[1,3,2]),[1,3,2]),3);

  [v, eigvs] = eig(T, CnM(:,:,1)); 
  [cc, idx] = max(abs(diag(eigvs)));
  x = v(:, idx);
  iter = iter+1;
  
end
tau = (1+lam)/(1-lam);

end