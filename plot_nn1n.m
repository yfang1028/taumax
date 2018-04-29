N = (fix(10.^[2:0.2:6]));

taupred = load('tauours_pred_nn1n');
taua = load('tauours_a_nn1n');
tauall = load('tauours_all_nn1n');
taupredacor = load('tauacor_pred_nn1n');
tauaacor = load('tauacor_a_nn1n');
tauallacor = load('tauacor_all_nn1n');
mserr = load('mserr_nn1n');
index = find(taupredacor<1);
taupredacor(index) = 1;

figure(1)
hold off
loglog(N, N./tauall,'k','linewidth',3.5)
hold on
loglog(N, N./taua, 'b--','linewidth',3.5)
loglog(N, N./taupred, 'r-.','linewidth',3.5)
loglog(N, 100*ones(1,length(N)),'color',[0.6,0.6,0.6],'linewidth',2.5)
set(gca, 'linewidth', 2, 'fontsize', 20)
xlabel('Number of Samples')
ylabel('Effective Sample Size')
legend('max','q1','prediction','100 ESS','location','north')
print -depsc nn1n_tauall.eps


figure(3)
hold off
loglog(N, tauall,'k','linewidth',3.5)
hold on
loglog(N, taua, 'b--','linewidth',3.5)
loglog(N, taupred, 'r-.', 'linewidth',3.5)
set(gca, 'linewidth', 2, 'fontsize', 20)
xlabel('Number of Samples')
ylabel('Autocorrelation Time')
legend('max','q1','prediction','100 ESS','location','northwest')
print -depsc nn1n_tauall_tau.eps




figure(5)
hold off
[sortx, sortidx] = sort(datax,'ascend');
plot(sortx,datay(sortidx),'k','marker','.','markersize',15,'linestyle','None')
hold on
tv1e2 = tt(1,sortidx);
plot(sortx,tv1e2,'k--','linewidth',3.5)
tv1e4 = tt(11,sortidx);
plot(sortx,tv1e4,'k-.','linewidth',3.5)
tv1e5 = tt(21,sortidx);
plot(sortx,tv1e5,'k','linewidth',3.5)
set(gca, 'linewidth', 2, 'fontsize', 20)
xlabel('x')
ylabel('y')
legend('data','N=1e2','N=1e4','N=1e6','location','southeast')
print -depsc nn1n_pred.eps

figure(6)
hold off
plot(ua(1:1e4),'k','linewidth',2.5)
axis([0,1e4,-4,4]);
set(gca, 'linewidth', 2, 'fontsize', 20)
xlabel('Sample Index')
ylabel('Function Value')
print -depsc nn1n_ua.eps

figure(7)
hold off
plot(uu(1:1e4),'k','linewidth',2.5)
axis([0,1e4,-4,4]);
set(gca, 'linewidth', 2, 'fontsize', 20)
xlabel('Sample Index')
ylabel('Function Value')
print -depsc nn1n_uall.eps

figure(8)
hold off
semilogx(N, mserr,'k','linewidth',2.5)
set(gca, 'linewidth', 2, 'fontsize', 20)
xlabel('Number of Samples')
ylabel('Mean Squared Error')
print -depsc nn1n_mse.eps


