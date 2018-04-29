warning ("off");
NMC = fix(10.^[2.0:0.5:6.0])-1;
s = length(NMC);
r = 12;
tauours_u1 = zeros(r, s);
tauours_u3 = zeros(r, s);
tauacor_u1 = zeros(r, s);
tauacor_u3 = zeros(r, s);
tauacor2_u1 = zeros(r, s);
tauacor2_u3 = zeros(r, s);
tauall_acor = zeros(r, s);
tauall_acor2 = zeros(r, s);
tauall_ours = zeros(r, s);

lam1 = exp(-0.1);
lam2 = exp(-0.2);
lam3 = exp(-0.3);
lam4 = exp(-0.4);
lam5 = exp(-0.5);
tau1 = (1+lam1)/(1-lam1);
tau2 = (1+lam2)/(1-lam2);
tau3 = (1+lam3)/(1-lam3);
tau4 = (1+lam4)/(1-lam4);
tau5 = (1+lam5)/(1-lam5);
taut = (tau1+tau2+tau3+tau4+tau5)/5;
taut_u1 = (tau1)/1;
taut_u2 = (tau1+tau2)/2;
taut_u3 = (tau1+tau2+tau3)/3;
taut_u4 = (tau1+tau2+tau3+tau4)/4;
taut_u5 = (tau1+tau2+tau3+tau4+tau5)/5;

c1 = lam1.^[1:max(NMC)-1]';
c2 = lam2.^[1:max(NMC)-1]';
c3 = lam3.^[1:max(NMC)-1]';
c4 = lam4.^[1:max(NMC)-1]';
c5 = lam5.^[1:max(NMC)-1]';
cu5 = (c1+c2+c3+c4+c5)/5;
cu4 = (c1+c2+c3+c4)/4;
cu3 = (c1+c2+c3)/3;
cu2 = (c1+c2)/2;
cu1 = (c1)/1;
cu = cu1;
ct = [1;cu];
posfix = '_u1';

for i=1:r
    tic;
    filename = ['1dsamples', num2str(i)];
    xs = load(filename);
    xs = xs(:, 2:end)';
    k = 1;
    for NN = NMC
        M = NN;
        N = NN;
        u1 = xs(:, 1:NN);
        phi5 = 1/120^0.5*(u1.^5-10*u1.^3+15*u1);
        phi4 = 1/24^0.5*(u1.^4-6*u1.^2+3);
        phi3 = 1/6^0.5*(u1.^3-3*u1);
        phi2 = 1/2^0.5*(u1.^2-1);
        phi1 = u1;
        u5 = phi1+phi2+phi3+phi4+phi5;
        u4 = phi1+phi2+phi3+phi4;
        u3 = phi1+phi2+phi3;
        u2 = phi1+phi2;
        u1 = phi1;
        
        % u1
        u = u1;
        [cf,uncf] = acffft(u, N);
        
        [lam, c0, sigma2, iter] = findlam2(uncf, 0, M-1);
        m1 = optm(lam, c0, sigma2);
        wours = lagwindowours(m1, M-1, lam);
        wours = wours';
        tauours_u1(i, k) = autocorfromw(cf', wours);
        [tauacor_u1(i, k), mmyacor] = myautocor(cf, M-1);
        [tauacor2_u1(i, k), mmyacor2] = myautocor2(cf, M-1);
        
        % u3
        u = u3;
        [cf, uncf] = acffft(u, N);
        
        [lam, c0, sigma2, iter] = findlam2(uncf, 0, M-1);
        m1 = optm(lam, c0, sigma2);
        wours = lagwindowours(m1, M-1, lam);
        wours = wours';
        tauours_u3(i, k) = autocorfromw(cf', wours);
        [tauacor_u3(i, k), mmyacor] = myautocor(cf, M-1);
        [tauacor2_u3(i, k), mmyacor] = myautocor2(cf, M-1);
        
        % worst case
        uall = [phi3+phi2+phi1; phi3-phi2+phi1; -phi3+phi2+phi1];
        CnM = cnmfft(uall, N);
        x0 = ones(3,1);

        [xacor, tauall_acor(i,k), cfx, m, w, iter] = optimizex6_acor(uall, CnM, x0, M-1);
        [xacor2, tauall_acor2(i,k), cfx, m, w, iter] = optimizex6_acor2(uall, CnM, x0, M-1);
        [xours, tau, cfx, uu, m, w, lam, c0, v, iter, isfail(i,k)] = optimizex6(uall, CnM, x0, M-1);
        if isfail(i,k) == 1
          tauall_ours(i,k) = tau;
        else
          [cf, uncf] = acffft(uu, N);        
          [lam, c0, sigma2, iter] = findlam2(uncf, 0, M-1);
          m1 = optm(lam, c0, sigma2);
          wours = lagwindowours(m1, M-1, lam);
          wours = wours';
          tauall_ours(i, k) = autocorfromw(cf', wours);
        end   
        k = k+1;
    end
    toc;
end

dlmwrite('1d_tauours_u1_new',tauours_u1)
dlmwrite('1d_tauours_u3_new',tauours_u3)
dlmwrite('1d_tauacor_u1_new',tauacor_u1)
dlmwrite('1d_tauacor2_u1_new',tauacor2_u1)
dlmwrite('1d_tauacor_u3_new',tauacor_u3)
dlmwrite('1d_tauacor2_u3_new',tauacor2_u3)
dlmwrite('1d_tauall_ours_new',tauall_ours)
dlmwrite('1d_tauall_acor_new',tauall_acor)
dlmwrite('1d_tauall_acor2_new',tauall_acor2)




