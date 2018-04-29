function w = lagwindowacor(cf, M)


w = zeros(1, M-1);
c = 10;
MM = fix(0.5*M/c);
m = 0;
for m = 2:MM
    tau = 1 + 2*sum(cf(2:m));
    if tau > 1 && m > tau*c
        w(1:m) = ones(1, m);
        break;
    end
end
if m == MM
    w(1:MM) = ones(1, MM);
end

