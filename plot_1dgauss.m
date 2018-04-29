taut = taut_u1;
plot_tauall_ours = load('1d_tauall_ours');
plot_tauall_acor = load('1d_tauall_ours');
N = 10.^[2:0.5:6];
mtauours = mean(plot_tauall_ours);
stauours = std(plot_tauall_ours);
mtauacor = mean(plot_tauall_acor);
stauacor = std(plot_tauall_acor);

figure(1)
hold off
semilogx(N, mtauours,'linewidth',2.5)
hold on
semilogx(N, mtauours-2.26*stauours/10^0.5,'r','linewidth',2.5)
semilogx(N, mtauours+2.26*stauours/10^0.5,'r','linewidth',2.5)
semilogx(N, ones(1,9)*taut,'k','linewidth',2.5)
axis([10^(2-0.05),10^(6),5,30])
set(gca, 'linewidth', 2, 'fontsize', 20)
xlabel('Number of Samples')
ylabel('Auto Correlation Time')
legend('Mean of Our Estimates','','95% confidence interval','true','location','northwest')
print -depsc 1dgauss_tauall.eps

taut = taut_u1;
plot_tauu1_ours = load('1d_tauours_u1');
plot_tauu1_acor = load('1d_tauacor_u1');
mtauours = mean(plot_tauu1_ours);
stauours = std(plot_tauu1_ours);
mtauacor = mean(plot_tauu1_acor);
stauacor = std(plot_tauu1_acor);
figure(2)
hold off
semilogx(N, mtauours,'k--','marker','o','linewidth',3.5)
hold on
semilogx(N, mtauacor,'k-.','marker','s','linewidth',3.5)
semilogx(N, stauours,'k--','linewidth',3.5)
semilogx(N, stauacor,'k-.','linewidth',3.5)
semilogx(N, ones(1,9)*taut,'k','linewidth',2.5)
axis([10^(2-0.05),10^6,0,30])
set(gca, 'linewidth', 2, 'fontsize', 20)
xlabel('Number of Samples')
ylabel('Mean/STD of Autocorrelation Time')
legend('new mean','acor mean','new STD','acor STD','true','location','northeast')
print -depsc 1dgauss_u1_compare.eps

taut = taut_u3;
plot_tauu3_ours = load('1d_tauours_u3');
plot_tauu3_acor = load('1d_tauacor_u3');
mtauours = mean(plot_tauu3_ours);
stauours = std(plot_tauu3_ours);
mtauacor = mean(plot_tauu3_acor);
stauacor = std(plot_tauu3_acor);
figure(3)
hold off

semilogx(N, mtauours,'k--','marker','o','linewidth',3.5)
hold on
semilogx(N, mtauacor,'k-.','marker','s','linewidth',3.5)
semilogx(N, stauours,'k--','linewidth',3.5)
semilogx(N, stauacor,'k-.','linewidth',3.5)
semilogx(N, ones(1,9)*taut,'k','linewidth',3.5)
axis([10^(2-0.05),10^6,0,20])
set(gca, 'linewidth', 2, 'fontsize', 20)
xlabel('Number of Samples')
ylabel('Mean/STD of Autocorrelation Time')
legend('new mean','acor mean','new STD','acor STD','true','location','northeast')
print -depsc 1dgauss_u3_compare.eps

taut = taut_u1;
plot_tauall_ours = load('1d_tauall_ours');
plot_tauall_acor = load('1d_tauall_acor');
mtauours = mean(plot_tauall_ours);
stauours = std(plot_tauall_ours);
mtauacor = mean(plot_tauall_acor);
stauacor = std(plot_tauall_acor);
figure(4)
hold off

semilogx(N, mtauours,'k--','marker','o','linewidth',3.5)
hold on
semilogx(N, mtauacor,'k-.','marker','s','linewidth',3.5)
semilogx(N, stauours,'k--','linewidth',3.5)
semilogx(N, stauacor,'k-.','linewidth',3.5)
plot(N, ones(1,9)*taut,'k','linewidth',3.5)
axis([10^(2-0.05), 10^6, 0, 30])
set(gca, 'linewidth', 2, 'fontsize', 20)
xlabel('Number of Samples')
ylabel('Mean/STD of Autocorrelation Time')
legend('new mean','acor mean','new STD','acor STD','true','location','northeast')
print -depsc 1dgauss_tauall_compare.eps


figure(5)
hold off
semilogx(N, ones(1,9)*taut,'k--','linewidth',2.5)
hold on
for i=1:12
    semilogx(N,plot_tauall_ours(i,:),'k','linewidth',2.5)
end

set(gca, 'linewidth', 2, 'fontsize', 20)
xlabel('Number of Samples')
ylabel('Autocorrelation Time')
legend('True value','Estimated value','location','northeast')
saveas(gca, '1dgauss_tauall_allrun.pdf')
print -depsc 1dgauss_tauall_allrun.eps
