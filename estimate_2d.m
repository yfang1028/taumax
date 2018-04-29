warning ("off");
NMC = fix(10.^[2:0.2:6])-1;
s = length(NMC);
r = 1;
tau1 = zeros(r, s);
tau2 = zeros(r, s);
tauall = zeros(r, s);
tau1acor = zeros(r, s);
tau2acor = zeros(r, s);
tauallacor = zeros(r, s);
isfail = zeros(r, s);
MINESS = 100;

for i=1:r
    tic;
    filename = ['2dsamples', num2str(i)];
    xs = load(filename);
    k = 1;
    for NN = NMC
        M = NN;
        N = NN;
        q1 = xs(1, 1:NN);
        q2 = xs(2, 1:NN);
        % q1
        [cf, uncf] = acffft(q1, N);        
        [lam, c0, sigma2, iter] = findlam2(uncf, 0, M-1);
        m1 = optm(lam, c0, sigma2);
        if m1<0 || imag(m1)~=0
            m1=0;
        end
        wours = lagwindowours(m1, M-1, lam);
        wours = wours';
        tau1(i, k) = autocorfromw(cf', wours);
        if (1+lam)/(1-lam) > N/MINESS
            tau1(i, k) = (1+lam)/(1-lam);
        end
        [tau1acor(i,k), mmyacor] = myautocor(cf, M-1);
        
        % q2
        [cf, uncf] = acffft(q2, N);        
        [lam, c0, sigma2, iter] = findlam2(uncf, 0, M-1);
        lamq2 = lam;
        m1 = optm(lam, c0, sigma2);
        if m1<0 || imag(m1)~=0
            m1=0;
        end       
        wours = lagwindowours(m1, M-1, lam);
        wours = wours';
        tau2(i, k) = autocorfromw(cf', wours);
        if (1+lam)/(1-lam) > N/MINESS
            tau2(i, k) = (1+lam)/(1-lam);
        end
        [tau2acor(i,k), mmyacor] = myautocor(cf, M-1);
        
        % worst
        qall = [q1; q2];
        CnM = cnmfft(qall, N);
        if tau1(i,k) > tau2(i,k)
            x0 = [1;0];
        else
            x0 = [0;1];
        end

        [xours, tau, cfx, uu, m, w, lam, c0, v, iter, isfail(i,k)] = optimizex6(qall, CnM, x0, M-1);
        if isfail(i,k) == 1
          tauall(i,k) = tau;
        else
          [cf,uncf] = acffft(uu, N);        
          [lam, c0, sigma2, iter] = findlam2(uncf, 0, M-1);
          m1 = optm(lam, c0, sigma2);
          if m1<0 || imag(m1)~=0
            m1=0;
          end
          wours = lagwindowours(m1, M-1, lam);
          wours = wours';
          tauall(i, k) = autocorfromw(cf', wours);
          if (1+lam)/(1-lam) > N/MINESS
                tauall(i, k) = (1+lam)/(1-lam);
            end
        end
        [xacor, tauallacor(i,k), cfx, macor, w, iter] = optimizex6_acor(qall, CnM, x0, M-1);        
        k = k+1;
    end
    toc;
end
dlmwrite('tau1_2d', tau1)
dlmwrite('tau2_2d', tau2)
dlmwrite('tauall_2d', tauall)
dlmwrite('tau1acor_2d', tau1acor)
dlmwrite('tau2acor_2d', tau2acor)
dlmwrite('tauallacor_2d', tauallacor)