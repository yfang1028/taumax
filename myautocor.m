function [tau, m] = myautocor(cf, M)

c = 10;
M = fix(0.5*M/c);

for m = 2:M
    tau = 1 + 2*sum(cf(2:m));
    if tau > 1 && m > tau*c
        break;
    end
end
if M < 2
    tau = 1;
else
    if m == M
        tau = 1 + 2*sum(cf(2:M-1));
    end
end










