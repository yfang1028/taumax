taupred = load('tauours_pred_logistic');
taua = load('tauours_a_logistic');
tauall = load('tauours_all_logistic');
taupredacor = load('tauacor_pred_logistic');
tauaacor = load('tauacor_a_logistic');
tauallacor = load('tauacor_all_logistic');
errortrain = load('errortrain_logistic');
errortest = load('errortest_logistic');
index = find(taupred<1);
taupred(index)=1;
index = find(taua<1);
taua(index)=1;
index = find(taupred<1.22223 & taupred>1.22221);
taupred(index)=1;
index = find(taua<1.22223 & taua>1.22221);
taua(index)=1;


N = (fix(10.^[1:0.2:5]));
figure(1)
hold off
loglog(N, N./tauall,'k-','linewidth',5.5)
hold on
loglog(N, N./taua, 'b--','linewidth',5.5)
loglog(N, N./taupred, 'r-.', 'linewidth',5.5)
loglog(N, 100*ones(1,length(N)),'color',[0.6,0.6,0.6],'linewidth',3.5)
set(gca, 'linewidth', 1, 'fontsize', 25,'fontname','times')
xlabel('Number of Samples')
ylabel('Effective Sample Size')
legend('Max','q_1','Pred','100 ESS','interpreter','latex','location','northwest')
saveas(gca, 'logistic_tauall.pdf')
print -depsc logistic_tauall.eps


figure(2)
hold off
loglog(N, tauall,'k','linewidth',5.5)
hold on
loglog(N, taua, 'b--','linewidth',5.5)
loglog(N, taupred, 'r-.', 'linewidth',5.5)
axis([10^1,10^5,10^(0-0.1),10^2])
set(gca, 'linewidth', 1, 'fontsize', 25,'fontname','times')
xlabel('Number of Samples')
ylabel('Autocorrelation Time')
legend('Max','q_1','Pred','interpreter','latex','location','northwest')
saveas(gca, 'logistic_tauall_tau.pdf')
print -depsc logistic_tauall_tau.eps

ndata = 690;
figure(3)
hold off
semilogx(N,errortrain/ndata*2,'k','linewidth',5.5)
hold on
semilogx(N,errortest/ndata*2,'r--','linewidth',5.5)
set(gca, 'linewidth', 1, 'fontsize', 25,'fontname','times')
xlabel('Number of Samples')
ylabel('Training/Testing Error')
legend('Training Error','Testing Error','location','northeast')
print -depsc logistic_error.eps





