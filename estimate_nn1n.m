warning("off");
load('datax');
load('datay');
NMC = fix(10.^[2:0.2:6.0])-1;
s = length(NMC);
r = 1;
taupred = zeros(r, s);
taua = zeros(r, s);
tauall = zeros(r, s);
taupredacor = zeros(r, s);
tauaacor = zeros(r, s);
tauallacor = zeros(r, s);
tauall_acor = zeros(r, s);
tauallacor_our = zeros(r, s);
isfail = zeros(r, s);
MINESS = 100;
mserr = zeros(r, s);
tt = zeros(s, 100);

beta = 2.5;
gamma = 0.8;

k = 1;
for i=1:r
    tic;
    filename = ['nn1nsamples',num2str(i)];
    xs = load(filename);
    k = 1;
    for NN = NMC
        M = NN;
        N = NN;
        uall = xs(:, 1:NN);
        ua = xs(1, 1:NN);
        [t, ts] = pred_nn1n(datax, uall);
        upred = ts(:, 1)';
        tt(k, :) = t;
        mserr(i, k) = sum((datay-t).^2)/100;
        % pred
        [cf, uncf] = acffft(upred, N);        
        [lam, c0, sigma2, iter] = findlam2(uncf, 0, M-1);
        m1 = optm(lam, c0, sigma2);
        if m1<0 || imag(m1)~=0
            m1 = 0;
        end
        wours = lagwindowours(m1, M-1, lam);
        wours = wours';
        taupred(i,k) = autocorfromw(cf', wours);
                    if (1+lam)/(1-lam) > N/MINESS
                taupred(i,k) = (1+lam)/(1-lam);
            end
        [taupredacor(i,k), mmyacor] = myautocor(cf, M-1);
        
        % a
        [cf, uncf] = acffft(ua, N);        
        [lam, c0, sigma2, iter] = findlam2(uncf, 0, M-1);
        m1 = optm(lam, c0, sigma2);
        if m1<0 || imag(m1)~=0
            m1 = 0;
        end
        wours = lagwindowours(m1, M-1, lam);
        wours = wours';
        taua(i, k) = autocorfromw(cf', wours);
        if (1+lam)/(1-lam) > N/MINESS
            taua(i, k) = (1+lam)/(1-lam);
        end
        [tauaacor(i, k), mmyacor] = myautocor(cf, M-1);
        
        % worst
        CnM = cnmfft(uall, N);
        for ui = 1:4
            [cf, uncf] = acffft(uall(ui,:), N);        
            [lam, c0, sigma2, iter] = findlam2(uncf, 0, M-1);
            m1 = optm(lam, c0, sigma2);
            if m1<0 || imag(m1)~=0
                m1 = 0;
            end
            wours = lagwindowours(m1, M-1, lam);
            wours = wours';
            tauui(ui) = autocorfromw(cf', wours);
            if (1+lam)/(1-lam) > N/MINESS
                tauui(ui) = (1+lam)/(1-lam);
            end
        end
        x0 = ones(4, 1);
        x0 = zeros(4, 1);
        [uia, uib] = max(tauui);
        x0(uib)=1;

        [xours, tau, cfx, uu, m, w, lam, c0, v, iter, isfail(i,k)] = optimizex6(uall, CnM, x0, M-1);
        if isfail(i,k) == 1
          tauall(i,k) = tau;
        else
          [cf, uncf] = acffft(uu, N);        
          [lam, c0, sigma2, iter] = findlam2(uncf, 0, M-1);
          m1 = optm(lam, c0, sigma2);
          if m1<0 || imag(m1)~=0
            m1=0;
          end
          wours = lagwindowours(m1, M-1, lam);
          wours = wours';
          tauall(i, k) = autocorfromw(cf', wours);
        end
        uu = xours'*uall;
        [cf, uncf] = acffft(uu, N);
        tauall_acor(i, k) = myautocor(cf, M-1);
        [xacor, tauallacor(i,k), cfx, macor, w, iter] = optimizex6_acor(uall, CnM, x0, M-1);
        uu = xacor'*uall;
        [cf, uncf] = acffft(uu, N); 
        tauallacor(i, k) = myautocor(cf, M-1);      
        [lam, c0, sigma2, iter] = findlam2(uncf, 0, M-1);
        m1 = optm(lam, c0, sigma2);
        wours = lagwindowours(m1, M-1, lam);
        wours = wours';
        tauallacor_our(i, k) = autocorfromw(cf', wours);
        
        k = k+1;
    end
    toc;
end
dlmwrite('tauours_pred_nn1n', taupred);
dlmwrite('tauours_a_nn1n', taua);
dlmwrite('tauours_all_nn1n', tauall);
dlmwrite('tauacor_pred_nn1n', taupredacor);
dlmwrite('tauacor_a_nn1n', tauaacor);
dlmwrite('tauacor_all_nn1n', tauallacor);
dlmwrite('tt_nn1n', tt);
dlmwrite('mserr_nn1n', mserr);
dlmwrite('tauours_all_acor_nn1n', tauall_acor);
dlmwrite('tauacor_all_our_nn1n', tauallacor_our);

