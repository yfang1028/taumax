function [lam, c0, sigma2, iter] = findlam2(c, m, M)
    lam = 1;
    lamprev = 0;
    iter = 0;
    while 1
        if abs(lam-lamprev)<1e-9 || iter>100 || lam>1
            break;
        end      
        k = m:M;  
        A = sum(c(k+1).*lam.^k);
        B = sum(k.*lam.^(2*k-1));
        C = sum(c(k+1).*k.*lam.^(k-1));
        D = sum(lam.^(2*k));
        Ap = sum(c(k+1).*k.*lam.^(k-1));
        Bp = sum(k.*(2*k-1).*lam.^(2*k-2));
        Cp = sum(c(k+1).*k.*(k-1).*lam.^(k-2));
        Dp = sum(2*k.*lam.^(2*k-1));       
        f = A*B-C*D;
        fp = A*Bp+Ap*B-C*Dp-D*Cp;
        lamprev = lam;
        lam = lam - f/fp;
        iter = iter+1;
    end
    c0 = A/D;
    ss = sum((c(k+1)-c0*lam.^k).^2);
    sigma2 = ss/(M-m+1);
    
    if lam > 1
        lam = lamprev;
    end
    
end