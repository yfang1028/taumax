function [tau, acf] = autocor(u, N, n)

isbias = 1;
[tau, acf] = acfs(u, N, n, isbias);
[m, M, lambda, c0, r] = window(acf);
w = lagwindow(m, M, lambda, c0, r, n);

tau = sum(acf.*w);
tau = tau + c0*lambda^(M+1)*(1-lambda^(n-M))/(1-lambda);

    