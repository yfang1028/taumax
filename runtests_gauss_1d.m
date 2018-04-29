%clear;
%rng('shuffle')
N = 1e6;
gap = 5; skip = 1;
nMC = gap*N; D = 1; dt = 1/50;
d = 1;
filename = 'temp';

x0 = zeros(d,1);
p0 = randn(d,1);

tic;

xs = BDsampler_gauss_1d(x0, dt, nMC, filename, gap, skip);

toc;
     
save 'xs_1d.mat' xs
