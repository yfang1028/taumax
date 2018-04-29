function CnM = cnmfft(u, N)

u = u(:, 1:N);
[d, N] = size(u);
for i = 1:d
  u(i, :) = u(i, :)- mean(u(i, :));
end
CnM = zeros(d, d, N);
f = zeros(d, 2*N);
for i = 1:d
    f(i, :) = fft(u(i, :), 2*N);
end
for i=1:d
    for j=1:d
        v(j, i, :) = real(ifft((f(i, :).*conj(f(j, :)))))/N;
    end
end
for i = 1:d
    for j = 1:d
        CnM(i, j, :) = 0.5*(v(i, j, 1:N)+v(j, i, 1:N));
    end
end 