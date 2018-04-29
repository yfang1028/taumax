function [cf,uncf] = acffft(u, N)

u = u(1:N);
f = fft(u - mean(u), 2*N);
uncf = real(ifft(f .* conj(f)))/N;
uncf = uncf(1:N);
cf = uncf/uncf(1);