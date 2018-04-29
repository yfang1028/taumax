tau1 = load('tau1_2d');
tau2 = load('tau2_2d');
tauall = load('tauall_2d');
tau1acor = load('tau1acor_2d');
tau2acor = load('tau2acor_2d');
tauallacor = load('tauallacor_2d');
N = (fix(10.^[2:0.1:6]));

figure(1)
hold off
loglog(N, tau1, 'k--', 'marker', 'o', 'markersize', 6, 'linewidth', 2.5)
hold on
loglog(N, tau2, 'k-.', 'marker', '>', 'markersize', 6, 'linewidth', 2.5)
axis([10^2, 10^6, 10^(0-0.1), 10^4])
set(gca, 'linewidth', 1, 'fontsize', 25)
xlabel('Number of Samples')
ylabel('Autocorrelation Time')
legend('q1', 'q2', '100 ESS', 'location', 'northwest')
print -depsc 2dgauss_tau12.eps

figure(2)
hold off
loglog(N, N./tau1, 'k--', 'marker', 'o', 'linewidth', 2.5)
hold on
loglog(N, N./tau2, 'k-.','marker','>', 'linewidth',2.5)
loglog(N, 100*ones(1,length(N)),'k','linewidth',2.5)
set(gca, 'linewidth', 2, 'fontsize', 20)
xlabel('Number of Samples')
ylabel('Effective Sample Size')
legend('q1','q2','100 ESS','location','northeast')
print -depsc 2dgauss_ess12.eps


figure(4)
hold off
scatter(q1(1:5000),q2(1:5000),10,'k',"filled")
set(gca, 'linewidth', 2, 'fontsize', 20)
xlabel('q1')
ylabel('q2')
axis([-2,6,-2,6])
print -depsc 2dgauss_scatter.eps

figure(5)
hold off
plot(q1(1:5000),'k','linewidth',3.5)
hold on
plot(q2(1:5000),'color',[0.6,0.6,0.6],'linewidth',3.5)
set(gca, 'linewidth', 2, 'fontsize', 20)
xlabel('Index of Samples')
ylabel('Function Values')
legend('q1','q2','location','northwest')
print -depsc 2dgauss_trajectory.eps

figure(6)
hold off
loglog(N, tauall,'k--','marker','s','linewidth',2.5)
hold on
loglog(N, tau1, 'k:','marker','o','linewidth',2.5)
loglog(N, tau2, 'k-.','marker','>', 'linewidth',2.5)
set(gca, 'linewidth', 2, 'fontsize', 20)
xlabel('Number of Samples')
ylabel('Autocorrelation Time')
legend('max','q1','q2','100 ESS','location','northwest')
saveas(gca, '2dgauss_tauall.pdf')
print -depsc 2dgauss_tauall.eps


figure(7)
hold off
loglog(N, N./tauall,'k--','marker','s','linewidth',2.5)
hold on
loglog(N, N./tau1, 'k:','marker','o','linewidth',2.5)
loglog(N, N./tau2, 'k-.','marker','>', 'linewidth',2.5)
loglog(N, 100*ones(1,length(N)),'k','linewidth',2.5)
set(gca, 'linewidth', 2, 'fontsize', 20)
xlabel('Number of Samples')
ylabel('Effective Sample Size')
legend('max','q1','q2','100 ESS','location','northwest')
saveas(gca, '2dgauss_essall.pdf')
print -depsc 2dgauss_essall.eps








