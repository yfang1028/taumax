function w = lagwindowours(m, M, lam)

w = zeros(1,M);
w(1:m) = ones(1, m);
w(m+1:M) = lam.^[1:M-m]';