warning("off");
datafile = 'australiandata';
data = load(datafile);
datax = data(:, 1:end-1);
ndata = size(datax, 1);
datastd = std(datax, 1);
datamean = mean(datax, 1);
datax = datax - repmat(datamean, ndata, 1);
datax = datax ./ repmat(datastd, ndata, 1);
datax = [datax, ones(ndata,1)];
datay = data(:, end);
datatrain = datax(1:2:end, :);
datatest = datax(2:2:end, :);
labeltrain = datay(1:2:end);
labeltest = datay(2:2:end);
NMC = fix(10.^[1:0.2:6.0])-1;
s = length(NMC);
r = 1;
taupred = zeros(r, s);
taua = zeros(r, s);
tauall = zeros(r, s);
taupredacor = zeros(r, s);
tauaacor = zeros(r, s);
tauallacor = zeros(r, s);
isfail = zeros(r, s);
errortrain = zeros(r, s);
errortest = zeros(r, s);
MINESS = 100;
d = size(datax, 2);

beta = 1;
gamma = 0.1;

k = 1;
for i=1:r
    tic;
    xs = load('xs_logistic');
    k = 1;
    for NN = NMC
        M = NN;
        N = NN;
        uall = xs(:, 1:NN);
        ua = xs(1, 1:NN);
        tstrain = datatrain * uall;
        tstrain = 1./(1+exp(-tstrain));
        ttrain = mean(tstrain, 2);
        tstest = datatest * uall;
        tstest = 1./(1+exp(-tstest));
        ttest = mean(tstest, 2);
        
        label = zeros(fix(ndata/2), 1);
        idx = find(ttrain>0.5);
        label(idx) = 1;
        errortrain(i, k) = sum(abs(label-labeltrain));
        
        label = zeros(fix(ndata/2), 1);
        idx = find(ttest>0.5);
        label(idx) = 1;
        errortest(i, k) = sum(abs(label-labeltest));
        upred = tstrain(1, :);
        % pred
        [cf, uncf] = acffft(upred, N);        
        [lam, c0, sigma2, iter] = findlam2(uncf, 0, M-1);
        if lam < 0
            lam = 0.1;
        end
        m1 = optm(lam, c0, sigma2);
        if m1 < 1
            m1 = 1;
        end
        wours = lagwindowours(m1, M-1, lam);
        wours = wours';
        taupred(i,k) = autocorfromw(cf', wours);
        if (1+lam)/(1-lam) > N/MINESS
            taupred(i, k) = (1+lam)/(1-lam);
        end
        [taupredacor(i, k), mmyacor] = myautocor(cf, M-1);
        
        % a
        [cf, uncf] = acffft(ua, N);        
        [lam, c0, sigma2, iter] = findlam2(uncf, 0, M-1);
        if lam<0
            lam=0.1;
        end
        m1 = optm(lam, c0, sigma2);
        if m1 < 1
            m1 = 1;
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
        MINESS = 100;
        tauui = zeros(1, d);
        for ui = 1:d
            [cf, uncf] = acffft(uall(ui, :), N);        
            [lam, c0, sigma2, iter] = findlam2(uncf, 0, M-1);
            if lam < 0
                lam = 0.1;
            end
            m1 = optm(lam, c0, sigma2);
            if m1 < 1
                m1 = 1;
            end
            wours = lagwindowours(m1, M-1, lam);
            wours = wours';
            tauui(ui) = autocorfromw(cf', wours);
            if (1+lam)/(1-lam) > N/MINESS
                tauui(ui) = (1+lam)/(1-lam);
            end
        end
        x0 = ones(d,1);
        x0 = zeros(d,1);
        [uia,uib] = max(tauui);
        x0(uib)=1;

        [xours, tau, cfx, uu, m, w, lam, c0, v, iter, isfail(i,k)] = optimizex6(uall, CnM, x0, M-1);
        if isfail(i,k) == 1
          tauall(i,k) = tau;
        else
          [cf, uncf] = acffft(uu, N);        
          [lam, c0, sigma2, iter] = findlam2(uncf, 0, M-1);
          m1 = optm(lam, c0, sigma2);
          wours = lagwindowours(m1, M-1, lam);
          wours = wours';
          tauall(i,k) = autocorfromw(cf', wours);
        end
        [xacor, tauallacor(i,k), cfx, macor, w, iter] = optimizex6_acor(uall, CnM, x0, M-1);
        
        k = k+1;
    end
    toc;
end
dlmwrite('tauours_pred_logistic', taupred);
dlmwrite('tauours_a_logistic', taua);
dlmwrite('tauours_all_logistic', tauall);
dlmwrite('tauacor_pred_logistic', taupredacor);
dlmwrite('tauacor_a_logistic', tauaacor);
dlmwrite('tauacor_all_logistic', tauallacor);
dlmwrite('errortrain_logistic', errortrain);
dlmwrite('errortest_logistic', errortest);
